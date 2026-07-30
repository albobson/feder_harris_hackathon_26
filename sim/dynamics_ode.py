"""Mean-field ODE model: non-lysogen (U_S, P_S phenotype split), lysogen (L), free
phage (P), driven by an exogenous piecewise-constant environment (sim/environment.py).

State vector y = [U_S, P_S, L, P_phage].

Integrated exactly segment-by-segment (constant coefficients within each environment
segment) rather than handing scipy.solve_ivp a single discontinuous RHS across the
whole trajectory -- this avoids the solver silently stepping over environment
switches. An optional logistic carrying-capacity term (-N/K, N = U_S+P_S+L) can be
added for bounded coexistence dynamics; without it (K=None) growth is unbounded,
which is fine for checking against the exact Lyapunov exponents in lib/fitness.py
(log(N(T))/T -> Lambda as T -> inf) but not for interior-equilibrium exploration.
"""
import numpy as np
from scipy.integrate import solve_ivp


def _rhs(t, y, state, params):
    U_S, P_S, L, P_phage = y
    p = params

    q = p["q_A"] if state == "A" else p["q_N"]
    g_S_U = p["g_A"] if state == "A" else p["g_N"] - p["D"]
    g_S_P = (p["g_A"] - p["c"]) if state == "A" else p["g_N"]
    g_L = p["g_A"] if state == "A" else p["g_N"]
    lam = p["lambda_A"] if state == "A" else p["lambda_N"]

    N = U_S + P_S + L
    density_term = -N / p["K"] if p["K"] is not None else 0.0

    switch_out_U = p["k_switch"] * q
    switch_out_P = p["k_switch"] * (1 - q)

    S_total = U_S + P_S
    infection = p["delta"] * P_phage * S_total
    # split infection loss proportionally between U_S and P_S
    frac_U = U_S / S_total if S_total > 0 else 0.0
    frac_P = P_S / S_total if S_total > 0 else 0.0

    dU_S = (g_S_U + density_term) * U_S - switch_out_U * U_S + switch_out_P * P_S - infection * frac_U
    dP_S = (g_S_P + density_term) * P_S - switch_out_P * P_S + switch_out_U * U_S - infection * frac_P
    dL = (g_L - lam + density_term) * L + infection
    dP_phage = p["beta"] * lam * L - p["delta"] * P_phage * S_total - p["m"] * P_phage

    return [dU_S, dP_S, dL, dP_phage]


def simulate_mean_field(segments, y0, k_switch, q_A, q_N, g_A, g_N, c, D,
                         lambda_A=0.0, lambda_N=0.0, delta=0.0, beta=0.0, m=1.0,
                         K=None, points_per_segment=5, rescale_threshold=1e50):
    """Integrate the mean-field ODE across all environment segments.

    Returns (t, y, log_scale) with y shape (4, len(t)): rows U_S, P_S, L, P_phage.
    `y` is periodically rescaled down whenever any component exceeds
    `rescale_threshold` (unbounded exponential growth otherwise overflows double
    precision long before a "long-run growth rate" trajectory is long enough to be
    informative). `log_scale[i]` is the cumulative log-factor removed by time
    `t[i]` -- add it back to `log(y)` to recover true (unscaled) log-population,
    e.g. for a long-run growth-rate estimate: `(log(y[-1,k]) + log_scale[-1] -
    log(y0[k])) / (t[-1]-t[0])`.
    """
    params = dict(k_switch=k_switch, q_A=q_A, q_N=q_N, g_A=g_A, g_N=g_N, c=c, D=D,
                  lambda_A=lambda_A, lambda_N=lambda_N, delta=delta, beta=beta, m=m, K=K)

    t_all = []
    y_all = []
    log_scale_all = []
    y_current = np.array(y0, dtype=float)
    cum_log_scale = 0.0

    for seg_start, seg_end, state in segments:
        if seg_end <= seg_start:
            continue
        t_eval = np.linspace(seg_start, seg_end, points_per_segment)
        sol = solve_ivp(_rhs, (seg_start, seg_end), y_current, args=(state, params),
                         t_eval=t_eval, method="LSODA", rtol=1e-8, atol=1e-10)
        if not sol.success:
            raise RuntimeError(f"ODE integration failed on segment {seg_start}-{seg_end}: {sol.message}")
        y_seg = sol.y
        if t_all and abs(t_eval[0] - t_all[-1][-1]) < 1e-12:
            t_all.append(t_eval[1:])
            y_all.append(y_seg[:, 1:])
        else:
            t_all.append(t_eval)
            y_all.append(y_seg)
        log_scale_all.append(np.full(y_all[-1].shape[1], cum_log_scale))
        y_current = y_seg[:, -1]

        max_val = np.max(np.abs(y_current))
        if max_val > rescale_threshold:
            factor = np.log(max_val)
            cum_log_scale += factor
            y_current = y_current / max_val

    t = np.concatenate(t_all)
    y = np.concatenate(y_all, axis=1)
    log_scale = np.concatenate(log_scale_all)
    return t, y, log_scale


def _sanity_checks():
    from environment import generate_environment_sequence

    rng = np.random.default_rng(1)
    # Fast-ish switching relative to t_max so the growth-rate estimate averages
    # over many environment cycles; t_max kept long enough to trigger the
    # rescaling path at least once (this deliberately tests the log_scale
    # bookkeeping, not just the no-rescale case).
    sigma_AN, sigma_NA = 0.2, 0.2
    t_max = 800.0
    segs = generate_environment_sequence(t_max, sigma_AN, sigma_NA, rng)

    # No phage (delta=beta=0): L should grow at exactly its Lyapunov exponent
    # Lambda_L (checked via lib/fitness.py) over a long trajectory, since L then
    # decouples from U_S/P_S entirely.
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "lib"))
    from fitness import lysogen_growth_rate

    g_A, g_N, lambda_A, lambda_N = 1.0, 0.4, 0.05, 0.05
    Lambda_L_exact = lysogen_growth_rate(sigma_AN, sigma_NA, g_A, g_N, lambda_A, lambda_N)

    y0 = [0.0, 0.0, 1.0, 0.0]  # start with 1 lysogen, no non-lysogens, no phage
    t, y, log_scale = simulate_mean_field(
        segs, y0, k_switch=1.0, q_A=0.1, q_N=0.95,
        g_A=g_A, g_N=g_N, c=0.3, D=5.0,
        lambda_A=lambda_A, lambda_N=lambda_N,
        delta=0.0, beta=0.0, m=1.0, K=None, rescale_threshold=1e6)
    assert log_scale[-1] > 0, "expected rescaling to have triggered at least once"
    L = y[2]
    true_log_L_final = np.log(L[-1]) + log_scale[-1]
    log_growth_rate = (true_log_L_final - np.log(L[0])) / (t[-1] - t[0])
    assert abs(log_growth_rate - Lambda_L_exact) < 0.02, (log_growth_rate, Lambda_L_exact)

    print("All sanity checks passed.")


if __name__ == "__main__":
    _sanity_checks()
