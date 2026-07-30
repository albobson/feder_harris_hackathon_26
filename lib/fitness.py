"""Growth-rate and invasion quantities for the torCAD bet-hedging / HK022-lysogeny
model. See derivation.md for the math and for two modeling errors that were caught
and corrected along the way.

Key structural point (this is what the model is *for*): lysogen and non-lysogen are
the SAME kind of object -- a population split between "unprepared" (U, torCAD off)
and "prepared" (P, torCAD on) cells, relaxing toward an environment-dependent target
prepared fraction q_E at rate k_switch. Growth rates depend only on phenotype and
environment, never on genotype. The genotypes differ in exactly two ways:

  1. their AEROBIC target prepared fraction q_A -- ~0.1 for the non-lysogen (the
     noisy bet-hedging state), 0 for the lysogen (uniformly off, which is precisely
     what Carey et al. 2019 Fig 5 says the prophage does);
  2. the lysogen additionally dies at the lysis rate lambda(E).

Their ANAEROBIC target q_N is the same, because the papers find lysogen and
non-lysogen anaerobic torCAD expression experimentally indistinguishable.

Modeling the lysogen this way -- rather than as a genotype whose expression tracks
oxygen instantaneously -- is what lets the model represent Carey Fig 1C/D, where it
is the LYSOGENS that fail to grow after rapid oxygen depletion (nobody is prepared
when the oxygen goes). An earlier version of this file gave the lysogen a lag-free,
penalty-free response, which silently handed it an advantage no cell has and
inverted the model's headline result. See derivation.md.
"""
import numpy as np
from scipy.linalg import expm


def stationary_fractions(sigma_AN, sigma_NA):
    """Stationary occupancy (p_A, p_N) of the two-state environment chain."""
    p_A = sigma_NA / (sigma_AN + sigma_NA)
    p_N = sigma_AN / (sigma_AN + sigma_NA)
    return p_A, p_N


def weighted_time_average(sigma_AN, sigma_NA, val_A, val_N):
    """p_A * val_A + p_N * val_N -- the mean rate for anything that is a simple
    (non-compartmental) function of the current environment, e.g. the lysis rate.
    """
    p_A, p_N = stationary_fractions(sigma_AN, sigma_NA)
    return p_A * val_A + p_N * val_N


def phenotype_growth_rates(state, g_A, g_N, c, D):
    """(unprepared, prepared) per-capita growth rates in the given environment.

    Depends on phenotype and environment ONLY -- not on genotype. Aerobically a
    prepared cell pays cost c for machinery it isn't using; anaerobically an
    unprepared cell pays the stall penalty D because it cannot respire TMAO yet.
    """
    if state == "A":
        return g_A, g_A - c
    return g_N - D, g_N


def _phenotype_matrix(state, k_switch, q_A, q_N, g_A, g_N, c, D, lam=0.0):
    """2x2 matrix for the (U, P) phenotype ODE dv/dt = A(state) v for one genotype,
    in the given (exogenous, currently-fixed) environment state.

    q_A / q_N are THIS genotype's target prepared fractions; lam is its lysis rate
    in this environment (0 for a non-lysogen), which applies uniformly to both
    phenotypes and so just shifts the whole matrix by -lam*I.
    """
    q = q_A if state == "A" else q_N
    g_U, g_P = phenotype_growth_rates(state, g_A, g_N, c, D)
    return np.array([
        [g_U - lam - k_switch * q, k_switch * (1 - q)],
        [k_switch * q, g_P - lam - k_switch * (1 - q)],
    ])


def _top_eigenvalue_2x2(M):
    tr = M[0, 0] + M[1, 1]
    det = M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]
    disc = tr * tr - 4 * det
    return (tr + np.sqrt(disc)) / 2.0


def growth_rate_numeric(segments, k_switch, q_A, q_N, g_A, g_N, c, D,
                         lambda_A=0.0, lambda_N=0.0, v0=(1.0, 1.0)):
    """Long-run growth rate of one genotype along a GIVEN realized environment path.

    This is the top Lyapunov exponent of a product of random matrices -- exact for
    any k_switch/sigma ratio, with no closed form in general (which is exactly why
    classical bet-hedging theory restricts its clean results to limiting regimes).
    Estimated directly: propagate a population vector through the matrix exponential
    of each environment segment, renormalizing every step and accumulating log-growth.

    Both genotypes MUST be evaluated on the SAME `segments` to be compared: they
    live in one shared environment, so the fair comparison is between their
    growth along one common realization, not between separate random draws.
    """
    v = np.array(v0, dtype=float)
    log_growth = 0.0
    t_total = 0.0
    max_rate = max(k_switch, abs(g_A), abs(g_A - c), abs(g_N), abs(g_N - D),
                   abs(lambda_A), abs(lambda_N), 1.0)
    chunk_dt = 200.0 / max_rate  # keeps any one exponent well under float64's ~709 limit

    for seg_start, seg_end, state in segments:
        dt_remaining = seg_end - seg_start
        if dt_remaining <= 0:
            continue
        lam = lambda_A if state == "A" else lambda_N
        A = _phenotype_matrix(state, k_switch, q_A, q_N, g_A, g_N, c, D, lam)
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


def growth_rate_fast_limit(sigma_AN, sigma_NA, k_switch, q_A, q_N, g_A, g_N, c, D,
                            lambda_A=0.0, lambda_N=0.0):
    """Closed form in the quasi-static-environment limit: each environment dwell is
    long compared to the phenotype system's own relaxation time, so within a dwell
    growth settles onto the TOP EIGENVALUE of the local frozen-environment matrix,
    and those combine across environments by a plain time-average.

    Note the validity condition is NOT simply k_switch >> sigma -- the local
    matrix's relaxation time depends on k_switch together with c, D, g_A, g_N. Using
    a naive q-weighted arithmetic mean of the bare compartment rates here (only
    correct as k_switch -> infinity) was an earlier bug; see derivation.md.
    """
    top_A = _top_eigenvalue_2x2(
        _phenotype_matrix("A", k_switch, q_A, q_N, g_A, g_N, c, D, lambda_A))
    top_N = _top_eigenvalue_2x2(
        _phenotype_matrix("N", k_switch, q_A, q_N, g_A, g_N, c, D, lambda_N))
    return weighted_time_average(sigma_AN, sigma_NA, top_A, top_N)


def growth_rate_slow_limit(sigma_AN, sigma_NA, g_A, g_N, c, D,
                            lambda_A=0.0, lambda_N=0.0):
    """Closed form in the k_switch -> 0 limit: phenotype is effectively a fixed
    lineage trait, so the population is dominated by whichever committed
    sub-lineage has the higher environment-averaged rate. This is the classical
    bet-hedging regime that most of the literature assumes -- and per NOTES.md it
    is probably NOT the relevant one here, since switching takes ~one generation.

    Independent of q_A/q_N: with no switching, the targets are never approached.
    """
    g_U_A, g_P_A = phenotype_growth_rates("A", g_A, g_N, c, D)
    g_U_N, g_P_N = phenotype_growth_rates("N", g_A, g_N, c, D)
    lam_bar = weighted_time_average(sigma_AN, sigma_NA, lambda_A, lambda_N)
    r_unprepared = weighted_time_average(sigma_AN, sigma_NA, g_U_A, g_U_N) - lam_bar
    r_prepared = weighted_time_average(sigma_AN, sigma_NA, g_P_A, g_P_N) - lam_bar
    return max(r_unprepared, r_prepared)


def phage_invasion_growth_rate(growth_advantage, lambda_bar, m, delta, S_star, beta):
    """Top eigenvalue of the linearized (rare lysogen + free phage) system against a
    resident non-lysogen population sitting at its ecological equilibrium S_star.

    `growth_advantage` is the lysogen's growth rate MINUS the resident non-lysogen's
    (both from growth_rate_numeric, both with lysis excluded -- lysis enters via
    lambda_bar). Measuring it relative to the resident is the whole point: in a
    density-regulated population the resident's own net per-capita growth is zero at
    equilibrium, so a rare lineage spreads iff it beats the RESIDENT, not iff it
    beats zero. An earlier version compared against zero, which put the invasion
    boundary roughly an order of magnitude too high in lambda. See derivation.md.

    Returns > 0 iff the phage-carrying lineage increases in frequency when rare.
    """
    A = np.array([
        [growth_advantage - lambda_bar, delta * S_star],
        [beta * lambda_bar, -(m + delta * S_star)],
    ])
    eigvals = np.linalg.eigvals(A)
    return float(np.max(eigvals.real))


def phage_invades(growth_advantage, lambda_bar, m, delta, S_star, beta):
    """Boolean form of the same condition, as a closed-form inequality:

        growth_advantage * (m + delta*S_star) > lambda_bar * (m - delta*S_star*(beta - 1))
    """
    lhs = growth_advantage * (m + delta * S_star)
    rhs = lambda_bar * (m - delta * S_star * (beta - 1))
    return lhs > rhs


def _sanity_checks():
    import sys
    import os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "sim"))
    from environment import generate_environment_sequence

    rng = np.random.default_rng(7)
    g_A, g_N, c, D = 1.0, 0.4, 0.3, 5.0
    q_A_NON, q_A_LYS, q_N = 0.1, 0.0, 0.95

    sigma = 0.2
    segs = generate_environment_sequence(20000.0, sigma, sigma, rng)

    # --- Limits agree with the numeric estimator in their regimes ---
    fast_pred = growth_rate_fast_limit(sigma, sigma, 1e4, q_A_NON, q_N, g_A, g_N, c, D)
    numeric_fast = growth_rate_numeric(segs, 1e4, q_A_NON, q_N, g_A, g_N, c, D)
    numeric_slow = growth_rate_numeric(segs, 1e-4, q_A_NON, q_N, g_A, g_N, c, D)
    slow_pred = growth_rate_slow_limit(sigma, sigma, g_A, g_N, c, D)
    assert abs(fast_pred - numeric_fast) < 0.02, (fast_pred, numeric_fast)
    assert abs(slow_pred - numeric_slow) < 0.02, (slow_pred, numeric_slow)
    assert abs(fast_pred - numeric_slow) > 0.005, "fast limit unexpectedly matches slow regime"

    # --- Lysis is a uniform tax: it shifts a genotype's growth rate by exactly
    # -lambda_bar, since it hits both phenotypes equally. ---
    lam = 0.07
    base = growth_rate_numeric(segs, 1.0, q_A_LYS, q_N, g_A, g_N, c, D)
    taxed = growth_rate_numeric(segs, 1.0, q_A_LYS, q_N, g_A, g_N, c, D,
                                 lambda_A=lam, lambda_N=lam)
    assert abs((base - taxed) - lam) < 1e-9, (base, taxed, lam)

    # --- The two genotypes differ ONLY via q_A, and with q_A equal they are
    # literally the same population. ---
    same_1 = growth_rate_numeric(segs, 1.0, 0.1, q_N, g_A, g_N, c, D)
    same_2 = growth_rate_numeric(segs, 1.0, 0.1, q_N, g_A, g_N, c, D)
    assert same_1 == same_2

    # --- The structural fix, stated as a test: with NO aerobic cost (c=0), being
    # uniformly off aerobically is never an advantage, because the only thing it
    # buys you is avoiding a cost that doesn't exist -- while it still leaves you
    # with nobody prepared when the oxygen goes. So the lysogen must lose (or at
    # best tie) at every environmental switching rate. The pre-fix model, which
    # let the lysogen respond instantly and skip the stall penalty D entirely,
    # got this backwards and had it WIN by a growing margin as sigma rose.
    for sig in [0.01, 0.1, 1.0, 10.0]:
        r = np.random.default_rng(3)
        sg = generate_environment_sequence(3000.0 / sig, sig, sig, r)
        lys = growth_rate_numeric(sg, 1.0, q_A_LYS, q_N, g_A, g_N, c=0.0, D=D)
        non = growth_rate_numeric(sg, 1.0, q_A_NON, q_N, g_A, g_N, c=0.0, D=D)
        assert lys <= non + 1e-9, f"sigma={sig}: lysogen {lys} should not beat non-lysogen {non} at c=0"

    # --- Phage invasion: with no horizontal transmission (delta=0) the phage's
    # fate collapses onto its host's relative advantage minus the lysis tax. ---
    rate = phage_invasion_growth_rate(0.2, 0.05, m=0.5, delta=0.0, S_star=2.0, beta=20.0)
    assert abs(rate - (0.2 - 0.05)) < 1e-10, rate

    # A lysogen with zero growth advantage and any lysis at all needs horizontal
    # transmission to survive; with beta=0 it strictly cannot invade.
    assert phage_invasion_growth_rate(0.0, 0.1, 0.5, 0.01, 50.0, 0.0) < 0

    # Boolean form matches the sign of the exact eigenvalue.
    rng2 = np.random.default_rng(0)
    for _ in range(2000):
        args = (rng2.uniform(-1, 1), rng2.uniform(0, 2), rng2.uniform(0.01, 2),
                rng2.uniform(0.01, 2), rng2.uniform(0.01, 5), rng2.uniform(0.1, 50))
        assert (phage_invasion_growth_rate(*args) > 0) == phage_invades(*args)

    print("All sanity checks passed.")


if __name__ == "__main__":
    _sanity_checks()
