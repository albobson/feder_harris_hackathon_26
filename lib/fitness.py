"""Closed-form and numerically-exact fitness quantities for the torCAD bet-hedging /
HK022-lysogeny model. See derivation.md for the math, what's exact vs. a limiting
approximation, and a corrected modeling error (environment is an exogenous, SHARED
condition -- not a per-cell compartment cells can be distributed across, unlike
phenotype, which genuinely is per-cell). All functions are plain and stateless.
"""
import numpy as np
from scipy.linalg import expm


def stationary_fractions(sigma_AN, sigma_NA):
    """Stationary occupancy (p_A, p_N) of the two-state environment chain."""
    p_A = sigma_NA / (sigma_AN + sigma_NA)
    p_N = sigma_AN / (sigma_AN + sigma_NA)
    return p_A, p_N


def weighted_time_average(sigma_AN, sigma_NA, val_A, val_N):
    """p_A * val_A + p_N * val_N -- the natural 'mean rate' for anything that is a
    simple (non-compartmental) function of the current environment.
    """
    p_A, p_N = stationary_fractions(sigma_AN, sigma_NA)
    return p_A * val_A + p_N * val_N


def lysogen_growth_rate(sigma_AN, sigma_NA, g_A, g_N, lambda_A=0.0, lambda_N=0.0):
    """Exact long-run growth rate of the lysogen population, for ANY switching
    rates. The lysogen's torCAD state is a deterministic function of the CURRENT
    environment (no internal phenotype compartment, no lag) -- so population size
    N obeys dN/dt = r(E(t))*N with r(A)=g_A-lambda_A, r(N)=g_N-lambda_N, and
    log(N(t)) is exactly the time-integral of r(E(s)). By the ergodic theorem this
    converges (almost surely, for any realized environment path) to the simple
    time-average, NOT to any eigenvalue of a "switching-rate matrix" -- there is no
    matrix here, because there is nothing to mix: environment is one shared
    condition, not a set of compartments the population is split across.
    """
    return weighted_time_average(sigma_AN, sigma_NA, g_A - lambda_A, g_N - lambda_N)


def _phenotype_matrix(state, k_switch, q_A, q_N, g_A, g_N, c, D):
    """2x2 matrix for the (U_S, P_S) phenotype-compartment ODE dv/dt = A(state) v,
    for the given (exogenous, currently-fixed) environment state.
    """
    if state == "A":
        q = q_A
        return np.array([
            [g_A - k_switch * q, k_switch * (1 - q)],
            [k_switch * q, (g_A - c) - k_switch * (1 - q)],
        ])
    else:
        q = q_N
        return np.array([
            [(g_N - D) - k_switch * q, k_switch * (1 - q)],
            [k_switch * q, g_N - k_switch * (1 - q)],
        ])


def _top_eigenvalue_2x2(M):
    tr = M[0, 0] + M[1, 1]
    det = M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]
    disc = tr * tr - 4 * det
    return (tr + np.sqrt(disc)) / 2.0


def nonlysogen_growth_rate_fast_limit(sigma_AN, sigma_NA, k_switch, q_A, q_N,
                                       g_A, g_N, c, D):
    """Closed-form growth rate in the "quasi-static environment" limit: each
    environment dwell is long compared to the (U_S,P_S) phenotype system's own
    relaxation time, so within a dwell the population's growth rate settles onto
    the TOP EIGENVALUE of the LOCAL (frozen-environment) 2x2 phenotype-mixing
    matrix -- then combine across environments with the same plain time-average
    used elsewhere (environment itself still isn't a compartment).

    An earlier version of this function used a naive q-weighted ARITHMETIC MEAN
    of the two phenotype compartments' bare growth rates instead of this local
    eigenvalue. That's only correct in the further k_switch -> infinity sub-limit
    (where the local eigenvalue itself converges to the naive average); at
    finite k_switch it was a poor approximation whenever k_switch wasn't ALSO
    large relative to the internal rate scales c, D, g_A, g_N (not just large
    relative to sigma) -- caught by comparing against
    nonlysogen_growth_rate_numeric on a case where D was large enough for the two
    formulas to disagree by >2x. This version is exact in the k_switch->infinity
    limit and a much better approximation than the old one at any finite
    k_switch, valid whenever environment dwell times are long relative to the
    LOCAL matrix's relaxation time (a combination of k_switch, c, D, g_A, g_N --
    not simply k_switch vs. sigma_AN/sigma_NA).
    """
    top_A = _top_eigenvalue_2x2(_phenotype_matrix("A", k_switch, q_A, q_N, g_A, g_N, c, D))
    top_N = _top_eigenvalue_2x2(_phenotype_matrix("N", k_switch, q_A, q_N, g_A, g_N, c, D))
    return weighted_time_average(sigma_AN, sigma_NA, top_A, top_N)


def nonlysogen_growth_rate_slow_limit(sigma_AN, sigma_NA, g_A, g_N, c, D):
    """Closed-form growth rate in the k_switch << sigma limit: phenotype is nearly
    a fixed lineage trait across many environment cycles (the classical bet-hedging
    regime most of the literature assumes -- see NOTES.md). The population is
    dominated by whichever committed sub-lineage (always-unprepared vs.
    always-prepared) has the higher environment-averaged growth rate; k_switch>0
    only matters for replenishing the losing type after a switch (not captured
    here -- this is the k_switch->0 growth-rate ceiling, not a full model of the
    replenishment dynamics).
    """
    r_unprepared = weighted_time_average(sigma_AN, sigma_NA, g_A, g_N - D)
    r_prepared = weighted_time_average(sigma_AN, sigma_NA, g_A - c, g_N)
    return max(r_unprepared, r_prepared)


def nonlysogen_growth_rate_numeric(segments, k_switch, q_A, q_N, g_A, g_N, c, D,
                                    v0=(1.0, 1.0)):
    """Numerically-exact top Lyapunov exponent of the (U_S, P_S) system for a
    GIVEN realized environment path (segments, from sim.environment), for any
    k_switch/sigma ratio. There is no closed form in the general case (a genuine
    'product of random matrices' problem -- the reason classical bet-hedging
    theory restricts itself to the fast/slow limits above). This estimates it
    directly: propagate a population vector through the exact matrix exponential
    of each environment segment, renormalizing every step to avoid overflow and
    accumulating the log-growth. Equivalent to (and a faster, non-stochastic
    alternative to) running sim/dynamics_ode.py with delta=beta=0 and reading off
    log(N(T))/T.
    """
    v = np.array(v0, dtype=float)
    log_growth = 0.0
    t_total = 0.0
    max_rate = max(k_switch, abs(g_A), abs(g_A - c), abs(g_N), abs(g_N - D), 1.0)
    chunk_dt = 200.0 / max_rate  # keeps any one exponent well under float64's ~709 limit

    for seg_start, seg_end, state in segments:
        dt_remaining = seg_end - seg_start
        if dt_remaining <= 0:
            continue
        A = _phenotype_matrix(state, k_switch, q_A, q_N, g_A, g_N, c, D)
        # Long dwell times need chunking -- a single expm(A*dt) with a large dt
        # overflows before we ever get a chance to renormalize.
        n_chunks = max(1, int(np.ceil(dt_remaining / chunk_dt)))
        dt_step = dt_remaining / n_chunks
        M_step = expm(A * dt_step)
        for _ in range(n_chunks):
            v = M_step @ v
            norm = np.linalg.norm(v)
            log_growth += np.log(norm)
            v = v / norm
        t_total += dt_remaining
    return log_growth / t_total


def phage_invasion_growth_rate(r_L, lambda_bar, m, delta, S_star, beta):
    """Top eigenvalue of the linearized (rare lysogen/phage) L-P system near an
    all-susceptible population at density S_star, using MEAN (time-averaged) host
    growth rate r_L and lysis rate lambda_bar -- i.e. this is itself a fast/
    quasi-static approximation in the environment, same spirit as the "fast
    limit" functions above (get r_L from lysogen_growth_rate(..., lambda_A=0,
    lambda_N=0) and lambda_bar from weighted_time_average(sigma_AN, sigma_NA,
    lambda_A, lambda_N)). Positive => the phage-carrying lineage grows when rare.
    """
    A = np.array([
        [r_L - lambda_bar, delta * S_star],
        [beta * lambda_bar, -(m + delta * S_star)],
    ])
    eigvals = np.linalg.eigvals(A)
    return float(np.max(eigvals.real))


def phage_invades(r_L, lambda_bar, m, delta, S_star, beta):
    """Boolean invasion condition, equivalent to phage_invasion_growth_rate(...) > 0
    but derived as a closed-form inequality (see derivation.md section 3):

        r_L * (m + delta*S_star)  >  lambda_bar * (m - delta*S_star*(beta - 1))
    """
    lhs = r_L * (m + delta * S_star)
    rhs = lambda_bar * (m - delta * S_star * (beta - 1))
    return lhs > rhs


def _sanity_checks():
    import sys
    import os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "sim"))
    from environment import generate_environment_sequence

    rng = np.random.default_rng(7)
    sigma_AN, sigma_NA = 0.2, 0.2
    g_A, g_N, lambda_A, lambda_N = 1.0, 0.4, 0.05, 0.05
    t_max = 20000.0
    segs = generate_environment_sequence(t_max, sigma_AN, sigma_NA, rng)

    # Lysogen: formula must match the direct path-integral of the instantaneous
    # rate to high precision (it's an exact ergodic average, not an approximation).
    Lambda_L = lysogen_growth_rate(sigma_AN, sigma_NA, g_A, g_N, lambda_A, lambda_N)
    log_N = sum((g_A - lambda_A if st == "A" else g_N - lambda_N) * (e - s) for s, e, st in segs)
    direct = log_N / t_max
    assert abs(Lambda_L - direct) < 0.01, (Lambda_L, direct)

    # Lysis strictly reduces lysogen growth rate.
    assert lysogen_growth_rate(sigma_AN, sigma_NA, g_A, g_N, 0.1, 0.1) < \
        lysogen_growth_rate(sigma_AN, sigma_NA, g_A, g_N, 0.0, 0.0)

    # Non-lysogen fast-limit closed form must match the numeric Lyapunov estimator
    # when environment dwell times are genuinely long relative to the LOCAL
    # matrix's own relaxation time (D=20 makes that relaxation time short, ~0.05,
    # so this needs sigma small enough that mean dwell 1/sigma is well above
    # that -- sigma=0.2 (mean dwell 5) is NOT long enough here and was exactly
    # the case that exposed the old naive-average formula being off by >2x).
    q_A, q_N, c, D = 0.1, 0.95, 0.3, 20.0
    slow_sigma = 0.001
    segs_slow_env = generate_environment_sequence(2_000_000.0, slow_sigma, slow_sigma, rng)
    numeric_at_k1 = nonlysogen_growth_rate_numeric(segs_slow_env, 1.0, q_A, q_N, g_A, g_N, c, D)
    fast_pred_at_k1 = nonlysogen_growth_rate_fast_limit(slow_sigma, slow_sigma, 1.0, q_A, q_N, g_A, g_N, c, D)
    assert abs(fast_pred_at_k1 - numeric_at_k1) < 0.02, (fast_pred_at_k1, numeric_at_k1)

    q_A, q_N, c, D = 0.1, 0.95, 0.3, 5.0
    fast_pred = nonlysogen_growth_rate_fast_limit(sigma_AN, sigma_NA, 1e4, q_A, q_N, g_A, g_N, c, D)
    numeric_fast = nonlysogen_growth_rate_numeric(segs, 1e4, q_A, q_N, g_A, g_N, c, D)
    numeric_slow = nonlysogen_growth_rate_numeric(segs, 1e-4, q_A, q_N, g_A, g_N, c, D)
    assert abs(fast_pred - numeric_fast) < 0.02, (fast_pred, numeric_fast)
    assert abs(fast_pred - numeric_slow) > 0.005, "fast limit unexpectedly matches slow regime too"

    # Slow-limit closed form should match the numeric estimator at tiny k_switch.
    slow_pred = nonlysogen_growth_rate_slow_limit(sigma_AN, sigma_NA, g_A, g_N, c, D)
    assert abs(slow_pred - numeric_slow) < 0.02, (slow_pred, numeric_slow)

    # Phage: no horizontal transmission (delta=0) collapses invasion rate onto the
    # host's own net rate (r_L - lambda_bar).
    r_L, lambda_bar, m, delta, S_star, beta = 0.2, 0.05, 0.5, 0.0, 2.0, 20.0
    rate = phage_invasion_growth_rate(r_L, lambda_bar, m, delta, S_star, beta)
    assert abs(rate - (r_L - lambda_bar)) < 1e-10, rate

    # phage_invades boolean matches sign of exact eigenvalue across a random sweep.
    rng2 = np.random.default_rng(0)
    for _ in range(2000):
        args = (rng2.uniform(0, 2), rng2.uniform(0, 2), rng2.uniform(0.01, 2),
                rng2.uniform(0.01, 2), rng2.uniform(0.01, 5), rng2.uniform(0.1, 50))
        assert (phage_invasion_growth_rate(*args) > 0) == phage_invades(*args)

    print("All sanity checks passed.")


if __name__ == "__main__":
    _sanity_checks()
