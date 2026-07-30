"""Mean-field ODE model of the two-genotype system, driven by an exogenous
piecewise-constant environment (sim/environment.py).

State y = [U_S, P_S, U_L, P_L, P_phage]: unprepared/prepared non-lysogens,
unprepared/prepared lysogens, free phage. Both genotypes carry a phenotype split and
relax at the same k_switch; they differ only in their aerobic target prepared
fraction (q_A_S vs q_A_L) and in that lysogens lyse. See lib/fitness.py for why that
symmetry matters -- an earlier version modelled the lysogen as a single compartment
tracking oxygen instantaneously, which let it skip the post-transition stall penalty
the non-lysogen paid.

Integrated exactly segment-by-segment (coefficients are constant within an
environment segment) rather than handing solve_ivp one discontinuous RHS, which
would let the solver step over environment switches. An optional logistic carrying
capacity K bounds total cell number; without it (K=None) growth is unbounded, which
is fine for reading off long-run growth rates but not for coexistence dynamics.
"""
import numpy as np
from scipy.integrate import solve_ivp

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "lib"))
from fitness import phenotype_growth_rates  # noqa: E402

U_S, P_S, U_L, P_L, PHAGE = 0, 1, 2, 3, 4


def _rhs(t, y, state, p):
    u_s, p_s, u_l, p_l, phage = y

    g_U, g_P = phenotype_growth_rates(state, p["g_A"], p["g_N"], p["c"], p["D"])
    q_S = p["q_A_S"] if state == "A" else p["q_N"]
    q_L = p["q_A_L"] if state == "A" else p["q_N"]
    lam = p["lambda_A"] if state == "A" else p["lambda_N"]
    k = p["k_switch"]

    n_cells = u_s + p_s + u_l + p_l
    dens = -n_cells / p["K"] if p["K"] is not None else 0.0

    inf_u = p["delta"] * phage * u_s
    inf_p = p["delta"] * phage * p_s

    du_s = (g_U + dens) * u_s - k * q_S * u_s + k * (1 - q_S) * p_s - inf_u
    dp_s = (g_P + dens) * p_s + k * q_S * u_s - k * (1 - q_S) * p_s - inf_p
    du_l = (g_U + dens - lam) * u_l - k * q_L * u_l + k * (1 - q_L) * p_l + inf_u
    dp_l = (g_P + dens - lam) * p_l + k * q_L * u_l - k * (1 - q_L) * p_l + inf_p
    dphage = p["beta"] * lam * (u_l + p_l) - inf_u - inf_p - p["m"] * phage

    return [du_s, dp_s, du_l, dp_l, dphage]


def simulate_mean_field(segments, y0, k_switch, q_A_S, q_A_L, q_N, g_A, g_N, c, D,
                         lambda_A=0.0, lambda_N=0.0, delta=0.0, beta=0.0, m=1.0,
                         K=None, points_per_segment=5, rescale_threshold=1e50):
    """Integrate across all environment segments.

    Returns (t, y, log_scale) with y shape (5, len(t)). `y` is rescaled down
    whenever any component exceeds `rescale_threshold` (unbounded exponential growth
    otherwise overflows float64 long before a trajectory is long enough to read a
    long-run growth rate off). `log_scale[i]` is the cumulative log-factor removed by
    time t[i]; add it back to log(y) to recover true log-population.
    """
    p = dict(k_switch=k_switch, q_A_S=q_A_S, q_A_L=q_A_L, q_N=q_N,
             g_A=g_A, g_N=g_N, c=c, D=D, lambda_A=lambda_A, lambda_N=lambda_N,
             delta=delta, beta=beta, m=m, K=K)

    t_all, y_all, scale_all = [], [], []
    y_cur = np.array(y0, dtype=float)
    cum_log_scale = 0.0

    for seg_start, seg_end, state in segments:
        if seg_end <= seg_start:
            continue
        t_eval = np.linspace(seg_start, seg_end, points_per_segment)
        sol = solve_ivp(_rhs, (seg_start, seg_end), y_cur, args=(state, p),
                         t_eval=t_eval, method="LSODA", rtol=1e-8, atol=1e-10)
        if not sol.success:
            raise RuntimeError(f"ODE failed on segment {seg_start}-{seg_end}: {sol.message}")
        if t_all and abs(t_eval[0] - t_all[-1][-1]) < 1e-12:
            t_all.append(t_eval[1:]); y_all.append(sol.y[:, 1:])
        else:
            t_all.append(t_eval); y_all.append(sol.y)
        scale_all.append(np.full(y_all[-1].shape[1], cum_log_scale))
        y_cur = sol.y[:, -1]

        peak = np.max(np.abs(y_cur))
        if peak > rescale_threshold:
            cum_log_scale += np.log(peak)
            y_cur = y_cur / peak

    return np.concatenate(t_all), np.concatenate(y_all, axis=1), np.concatenate(scale_all)


def _sanity_checks():
    sys.path.insert(0, os.path.dirname(__file__))
    from environment import generate_environment_sequence
    from fitness import growth_rate_numeric

    rng = np.random.default_rng(1)
    sigma = 0.2
    segs = generate_environment_sequence(600.0, sigma, sigma, rng)
    g_A, g_N, c, D, k = 1.0, 0.4, 0.1, 5.0, 1.0
    q_A_S, q_A_L, q_N = 0.1, 0.0, 0.95

    # Cross-validate the ODE layer against the analytic layer: with no phage
    # transmission each genotype grows at its own Lyapunov exponent, and the two
    # implementations are completely independent routes to the same number.
    for label, idx_pair, q_A, lam in [
        ("non-lysogen", (U_S, P_S), q_A_S, 0.0),
        ("lysogen", (U_L, P_L), q_A_L, 0.06),
    ]:
        y0 = np.zeros(5)
        y0[idx_pair[0]] = 1.0
        t, y, log_scale = simulate_mean_field(
            segs, y0, k_switch=k, q_A_S=q_A_S, q_A_L=q_A_L, q_N=q_N,
            g_A=g_A, g_N=g_N, c=c, D=D, lambda_A=lam, lambda_N=lam,
            delta=0.0, beta=0.0, m=1.0, K=None, rescale_threshold=1e6)
        tot = y[idx_pair[0]] + y[idx_pair[1]]
        ode_rate = (np.log(tot[-1]) + log_scale[-1] - np.log(tot[0])) / (t[-1] - t[0])
        analytic = growth_rate_numeric(segs, k, q_A, q_N, g_A, g_N, c, D,
                                        lambda_A=lam, lambda_N=lam)
        assert abs(ode_rate - analytic) < 0.02, (label, ode_rate, analytic)

    # The structural fix, visible in the ODE layer too: with c=0 the lysogen's
    # only distinguishing feature is having nobody prepared aerobically, which is
    # a pure liability, so it must not out-grow the non-lysogen.
    y0 = np.array([1.0, 0.0, 1.0, 0.0, 0.0])
    t, y, _ = simulate_mean_field(
        segs, y0, k_switch=k, q_A_S=q_A_S, q_A_L=q_A_L, q_N=q_N,
        g_A=g_A, g_N=g_N, c=0.0, D=D, delta=0.0, beta=0.0, m=1.0,
        K=None, rescale_threshold=1e6)
    assert (y[U_L] + y[P_L])[-1] <= (y[U_S] + y[P_S])[-1] * (1 + 1e-9)

    # Carrying capacity bounds the population.
    y0 = np.array([10.0, 0.0, 0.0, 0.0, 0.0])
    _, y, _ = simulate_mean_field(
        segs, y0, k_switch=k, q_A_S=q_A_S, q_A_L=q_A_L, q_N=q_N,
        g_A=g_A, g_N=g_N, c=c, D=D, delta=0.0, beta=0.0, m=1.0, K=1000.0)
    assert np.max(y[U_S] + y[P_S]) < 5000.0

    print("All sanity checks passed.")


if __name__ == "__main__":
    _sanity_checks()
