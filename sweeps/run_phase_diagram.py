"""Parameter-grid driver for the host/phage phase diagram.

Sweeps (sigma, k_switch, c, lambda0, alpha) and computes, at every grid point, BOTH
a host outcome (does the deterministic-switching lysogen out-grow the bet-hedging
non-lysogen?) and a phage outcome (does the phage invade when rare?) -- reported as
separate columns, never combined into one score (see NOTES.md).

Uses the fast analytical/numeric layer (lib/fitness.py), not the stochastic
simulator -- appropriate for dense grids; the stochastic simulator
(sim/dynamics_stochastic.py, via sim/metrics.py) is for specific finite-population
questions (does a SPECIFIC rare mutant establish), not for scanning thousands of
grid points.
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "lib"))
from fitness import (  # noqa: E402
    lysogen_growth_rate,
    nonlysogen_growth_rate_fast_limit,
    nonlysogen_growth_rate_numeric,
    nonlysogen_growth_rate_slow_limit,
    phage_invasion_growth_rate,
    weighted_time_average,
)

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "sim"))
from environment import generate_environment_sequence  # noqa: E402


def _env_cache_key(sigma_AN, sigma_NA, seed):
    return (round(sigma_AN, 12), round(sigma_NA, 12), seed)


def _compute_point(sigma_AN, sigma_NA, k_switch, q_A, q_N, g_A, g_N, c, D,
                    lambda0, alpha, m, delta, beta, S_star, segments):
    lambda_A = lambda0 * (1 + alpha)
    lambda_N = lambda0

    Lambda_L = lysogen_growth_rate(sigma_AN, sigma_NA, g_A, g_N, lambda_A, lambda_N)
    Lambda_S = nonlysogen_growth_rate_numeric(segments, k_switch, q_A, q_N, g_A, g_N, c, D)
    host_gap = Lambda_L - Lambda_S  # >0: lysogeny wins the host-growth comparison

    r_L = lysogen_growth_rate(sigma_AN, sigma_NA, g_A, g_N, 0.0, 0.0)
    lambda_bar = weighted_time_average(sigma_AN, sigma_NA, lambda_A, lambda_N)
    phage_rate = phage_invasion_growth_rate(r_L, lambda_bar, m, delta, S_star, beta)

    return dict(sigma_AN=sigma_AN, sigma_NA=sigma_NA, k_switch=k_switch, c=c, D=D,
                lambda0=lambda0, alpha=alpha,
                Lambda_L=Lambda_L, Lambda_S=Lambda_S, host_gap=host_gap,
                phage_invasion_rate=phage_rate)


def run_phase_diagram(sigma_values, k_switch_values, c_values, lambda0_values,
                       alpha_values, q_A, q_N, g_A, g_N, D, m, delta, beta, S_star,
                       sigma_NA_values=None, target_env_segments=3000, env_seed=0,
                       n_jobs=-1):
    """Returns a pandas DataFrame, one row per grid point.

    sigma_values is used for BOTH sigma_AN and sigma_NA unless sigma_NA_values is
    given explicitly (same length as sigma_values, paired index-wise) -- most
    exploration will want a single "environmental switching rate" axis with
    symmetric dwell times, but asymmetric dwell (e.g. long anaerobic, short
    aerobic bouts) is supported via sigma_NA_values.

    Each sigma value gets its own environment realization sized so it has
    roughly `target_env_segments` segments regardless of sigma (mean dwell time
    is 1/sigma, so total simulated time = target_env_segments / sigma) -- a
    single FIXED total time across a wide sigma range is a trap: at low sigma
    it under-samples (too few segments to converge), and at high sigma it
    over-samples catastrophically (a sigma=10 point with a fixed t_max=200,000
    generates ~2 million segments, both slow to compute and expensive to ship
    to worker processes for every grid point sharing that sigma).
    """
    if sigma_NA_values is None:
        sigma_NA_values = sigma_values

    # One shared environment realization per sigma pair, reused across all other
    # (k_switch, c, lambda0, alpha) grid points at that sigma -- environment
    # doesn't depend on those, and generating a fresh path per point would be
    # both wasteful and would inject spurious point-to-point noise into what
    # should be a smooth sweep.
    env_by_sigma = {}
    for sigma_AN, sigma_NA in zip(sigma_values, sigma_NA_values):
        rng = np.random.default_rng(env_seed)
        env_t_max = target_env_segments / max(sigma_AN, sigma_NA)
        env_by_sigma[_env_cache_key(sigma_AN, sigma_NA, env_seed)] = \
            generate_environment_sequence(env_t_max, sigma_AN, sigma_NA, rng)

    # Build the flat argument list with each point's segments passed explicitly
    # (not looked up from a shared dict captured by closure) -- keeps what gets
    # pickled to each worker down to the one relevant segment list.
    jobs = []
    for sigma_AN, sigma_NA in zip(sigma_values, sigma_NA_values):
        segments = env_by_sigma[_env_cache_key(sigma_AN, sigma_NA, env_seed)]
        for k_switch, c, lambda0, alpha in itertools.product(
                k_switch_values, c_values, lambda0_values, alpha_values):
            jobs.append((sigma_AN, sigma_NA, k_switch, c, lambda0, alpha, segments))

    def _job(item):
        sigma_AN, sigma_NA, k_switch, c, lambda0, alpha, segments = item
        return _compute_point(sigma_AN, sigma_NA, k_switch, q_A, q_N, g_A, g_N, c, D,
                               lambda0, alpha, m, delta, beta, S_star, segments)

    results = Parallel(n_jobs=n_jobs)(delayed(_job)(item) for item in jobs)
    return pd.DataFrame(results)


def _sanity_checks():
    # NOTE on scope: this checks the sweep's PLUMBING (grid shape, and that its
    # per-point calculation matches calling lib.fitness directly with the same
    # inputs), plus the one qualitative pattern that's already been independently
    # verified in lib/fitness.py (fast vs. slow k_switch limits, at FIXED sigma).
    # It deliberately does NOT assert a specific sigma-driven (in-host vs.
    # free-living) direction for host_gap: getting the fast/slow k_switch limits
    # to disagree in sign required deliberately hand-tuned parameters (see
    # derivation.md's fast-limit correction) and isn't something to bake into a
    # plumbing test as if it were settled -- that qualitative question is the
    # open scientific question this model exists to explore, not a known answer.
    sigma_AN, sigma_NA = 0.05, 0.05
    q_A, q_N, g_A, g_N, c, D = 0.1, 0.95, 1.0, 0.4, 0.3, 5.0
    lambda0, alpha, m, delta, beta, S_star = 0.02, 0.0, 0.5, 0.01, 20.0, 50.0

    df = run_phase_diagram(
        sigma_values=[sigma_AN], k_switch_values=[1e-4, 1.0, 1e4],
        c_values=[c], lambda0_values=[lambda0], alpha_values=[alpha],
        q_A=q_A, q_N=q_N, g_A=g_A, g_N=g_N, D=D, m=m, delta=delta, beta=beta,
        S_star=S_star, target_env_segments=2000, env_seed=0, n_jobs=1,
    )
    assert len(df) == 3
    assert set(df.columns) == {"sigma_AN", "sigma_NA", "k_switch", "c", "D",
                                "lambda0", "alpha", "Lambda_L", "Lambda_S",
                                "host_gap", "phage_invasion_rate"}

    # Plumbing check: recompute Lambda_L directly and confirm it matches every row
    # (Lambda_L doesn't depend on k_switch, so all three rows should agree).
    lambda_A_expected = lambda0 * (1 + alpha)
    Lambda_L_direct = lysogen_growth_rate(sigma_AN, sigma_NA, g_A, g_N, lambda_A_expected, lambda0)
    assert np.allclose(df["Lambda_L"].values, Lambda_L_direct, atol=1e-9)

    # Already-verified qualitative pattern (see lib/fitness.py, derivation.md):
    # at fixed sigma, Lambda_S at very small k_switch should sit near the
    # slow-limit closed form, and at very large k_switch near the (corrected,
    # local-eigenvalue) fast-limit closed form -- these two limits need not be
    # equal, but Lambda_S should land close to each in the corresponding regime.
    df_by_k = df.set_index("k_switch")
    slow_pred = nonlysogen_growth_rate_slow_limit(sigma_AN, sigma_NA, g_A, g_N, c, D)
    fast_pred = nonlysogen_growth_rate_fast_limit(sigma_AN, sigma_NA, 1e4, q_A, q_N, g_A, g_N, c, D)
    assert abs(df_by_k.loc[1e-4, "Lambda_S"] - slow_pred) < 0.02
    assert abs(df_by_k.loc[1e4, "Lambda_S"] - fast_pred) < 0.02

    print(df[["k_switch", "Lambda_L", "Lambda_S", "host_gap", "phage_invasion_rate"]])
    print("All sanity checks passed.")


if __name__ == "__main__":
    _sanity_checks()
