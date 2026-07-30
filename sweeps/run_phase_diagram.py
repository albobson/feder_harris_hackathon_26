"""Parameter-grid driver for the host/phage phase diagram.

At every grid point computes BOTH outcomes, as separate columns never combined into
one score (see NOTES.md):

  host_gap             Lambda_L - Lambda_S, the lysogen's growth-rate advantage over
                       the bet-hedging non-lysogen. Both genotypes are evaluated with
                       the SAME phenotype machinery on the SAME realized environment
                       path, differing only in aerobic target prepared fraction
                       (q_A_L vs q_A_S) and in the lysogen's lysis rate.
  phage_invasion_rate  whether the phage-carrying lineage increases in FREQUENCY when
                       rare, i.e. measured against the resident non-lysogen's growth
                       rate, not against zero.

Uses the analytic/numeric layer (lib/fitness.py). The stochastic simulator
(sim/dynamics_stochastic.py via sim/metrics.py) answers finite-population questions
about specific points, and is far too slow for dense grids.
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "lib"))
from fitness import (  # noqa: E402
    growth_rate_numeric,
    phage_invasion_growth_rate,
    weighted_time_average,
)

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "sim"))
from environment import generate_environment_sequence  # noqa: E402


def _env_key(sigma_AN, sigma_NA, seed):
    return (round(sigma_AN, 12), round(sigma_NA, 12), seed)


def _compute_point(sigma_AN, sigma_NA, k_switch, q_A_S, q_A_L, q_N, g_A, g_N, c, D,
                    lambda0, alpha, m, delta, beta, S_star, segments):
    lambda_A = lambda0 * (1 + alpha)
    lambda_N = lambda0

    # Both genotypes on the same environment path, same machinery, same k_switch.
    Lambda_S = growth_rate_numeric(segments, k_switch, q_A_S, q_N, g_A, g_N, c, D)
    Lambda_L_nolysis = growth_rate_numeric(segments, k_switch, q_A_L, q_N, g_A, g_N, c, D)
    Lambda_L = growth_rate_numeric(segments, k_switch, q_A_L, q_N, g_A, g_N, c, D,
                                    lambda_A=lambda_A, lambda_N=lambda_N)
    host_gap = Lambda_L - Lambda_S  # > 0: the on/off switching genotype out-grows

    # Phage invasion is judged relative to the resident, and the lysis tax is
    # carried by lambda_bar rather than folded into the growth advantage.
    growth_advantage = Lambda_L_nolysis - Lambda_S
    lambda_bar = weighted_time_average(sigma_AN, sigma_NA, lambda_A, lambda_N)
    phage_rate = phage_invasion_growth_rate(growth_advantage, lambda_bar, m, delta,
                                             S_star, beta)

    return dict(sigma_AN=sigma_AN, sigma_NA=sigma_NA, k_switch=k_switch, c=c, D=D,
                lambda0=lambda0, alpha=alpha,
                Lambda_S=Lambda_S, Lambda_L=Lambda_L,
                growth_advantage=growth_advantage, host_gap=host_gap,
                phage_invasion_rate=phage_rate)


def run_phase_diagram(sigma_values, k_switch_values, c_values, lambda0_values,
                       alpha_values, q_A_S, q_A_L, q_N, g_A, g_N, D, m, delta, beta,
                       S_star, sigma_NA_values=None, target_env_segments=3000,
                       env_seed=0, n_jobs=-1):
    """Returns a DataFrame, one row per grid point.

    sigma_values is used for both sigma_AN and sigma_NA unless sigma_NA_values is
    given explicitly (paired index-wise) -- most exploration wants one symmetric
    "environmental switching rate" axis, but asymmetric dwell times are supported.

    Each sigma gets an environment realization sized to hold roughly
    `target_env_segments` segments (total time = target/sigma). A single fixed total
    time across a wide sigma range is a trap: it under-samples at low sigma and
    over-samples catastrophically at high sigma.
    """
    if sigma_NA_values is None:
        sigma_NA_values = sigma_values

    env_by_sigma = {}
    for sigma_AN, sigma_NA in zip(sigma_values, sigma_NA_values):
        rng = np.random.default_rng(env_seed)
        t_max = target_env_segments / max(sigma_AN, sigma_NA)
        env_by_sigma[_env_key(sigma_AN, sigma_NA, env_seed)] = \
            generate_environment_sequence(t_max, sigma_AN, sigma_NA, rng)

    jobs = []
    for sigma_AN, sigma_NA in zip(sigma_values, sigma_NA_values):
        segments = env_by_sigma[_env_key(sigma_AN, sigma_NA, env_seed)]
        for k_switch, c, lambda0, alpha in itertools.product(
                k_switch_values, c_values, lambda0_values, alpha_values):
            jobs.append((sigma_AN, sigma_NA, k_switch, c, lambda0, alpha, segments))

    def _job(item):
        sigma_AN, sigma_NA, k_switch, c, lambda0, alpha, segments = item
        return _compute_point(sigma_AN, sigma_NA, k_switch, q_A_S, q_A_L, q_N,
                               g_A, g_N, c, D, lambda0, alpha, m, delta, beta,
                               S_star, segments)

    results = Parallel(n_jobs=n_jobs)(delayed(_job)(item) for item in jobs)
    return pd.DataFrame(results)


def _sanity_checks():
    q_A_S, q_A_L, q_N = 0.1, 0.0, 0.95
    g_A, g_N, c, D = 1.0, 0.4, 0.1, 5.0

    df = run_phase_diagram(
        sigma_values=[0.01, 0.1, 1.0, 10.0], k_switch_values=[1.0], c_values=[c],
        lambda0_values=[0.0, 0.2], alpha_values=[0.0],
        q_A_S=q_A_S, q_A_L=q_A_L, q_N=q_N, g_A=g_A, g_N=g_N, D=D,
        m=1.0, delta=0.001, beta=3.0, S_star=10.0,
        target_env_segments=3000, env_seed=1, n_jobs=-1)
    assert len(df) == 8

    # Lysis is a pure tax on the lysogen: host_gap must fall by exactly lambda_bar.
    for sigma in [0.01, 0.1, 1.0, 10.0]:
        sub = df[df["sigma_AN"] == sigma].set_index("lambda0")
        assert abs((sub.loc[0.0, "host_gap"] - sub.loc[0.2, "host_gap"]) - 0.2) < 1e-9

    # The corrected structural result: with the lysogen no longer given a lag-free
    # response, its advantage DECREASES as the environment gets more volatile, and
    # goes negative. The pre-fix model had this backwards (advantage growing with
    # volatility) and SUMMARY.md reported that artifact as a real finding.
    #
    # The decline saturates rather than continuing forever: once sigma is well past
    # k_switch neither genotype can track the environment at all, both settle to a
    # fixed composition, and host_gap flattens out (verified on a finer grid --
    # it bottoms out near sigma ~ 3 and stays there). So assert the sign change and
    # a substantial decline, and allow the plateau, rather than strict monotonicity.
    gaps = df[df["lambda0"] == 0.0].sort_values("sigma_AN")["host_gap"].values
    assert gaps[0] > 0 > gaps[-1], f"expected a sign change across sigma, got {gaps}"
    assert gaps[-1] < gaps[0] - 0.02, f"expected a substantial decline, got {gaps}"
    assert np.all(np.diff(gaps) < 5e-3), f"expected no real rebound in sigma, got {gaps}"

    print(df[["sigma_AN", "lambda0", "Lambda_S", "Lambda_L", "host_gap",
              "phage_invasion_rate"]].to_string(index=False))
    print("All sanity checks passed.")


if __name__ == "__main__":
    _sanity_checks()
