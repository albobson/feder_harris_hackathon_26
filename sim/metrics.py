"""Outcome metrics computed from the stochastic simulator, plus a thin decomposition
helper over lib/fitness.py.

Host and phage outcomes stay SEPARATE metrics throughout (see NOTES.md) -- a caller
assembling a phase diagram should plot both, not a combined score.

Everything here answers finite-population questions the deterministic layer cannot:
whether a lineage founded by a single cell survives its first few divisions, and
whether a phage introduction persists. Those depend on birth and death rates
separately, so they are sensitive to the baseline turnover `d0` (see
dynamics_stochastic.py) in a way that growth-rate comparisons are not.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "lib"))
from fitness import phage_invasion_growth_rate  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
from environment import generate_environment_sequence  # noqa: E402
from dynamics_stochastic import (  # noqa: E402
    simulate_stochastic, U_S, P_S, U_L, P_L, PHAGE,
)


def _run_replicates(rng, n_reps, t_max, sigma_AN, sigma_NA, y0, success, failure,
                     **sim_kwargs):
    """Fraction of replicates hitting `success(counts)` before `failure(counts)`.

    Both predicates are checked as stopping conditions, so an established lineage
    never has to be simulated past the point where the outcome is settled.
    """
    n_success = 0
    for _ in range(n_reps):
        child = np.random.default_rng(rng.integers(0, 2**63 - 1))
        segs = generate_environment_sequence(t_max, sigma_AN, sigma_NA, child)
        _, hist = simulate_stochastic(
            segs, y0, rng=child,
            extinction_check=lambda ct: success(ct) or failure(ct), **sim_kwargs)
        n_success += bool(success(hist[-1]))
    return n_success / n_reps


def lysogen_establishment_probability(rng, sigma_AN, sigma_NA, k_switch, q_A_S,
                                       q_A_L, q_N, g_A, g_N, c, D,
                                       lambda_A=0.0, lambda_N=0.0, d0=0.5,
                                       threshold=25, n_reps=300, t_max=200.0):
    """P(a lineage founded by ONE lysogen cell reaches `threshold` cells).

    No horizontal transmission: this isolates whether the lysogen genotype's own
    demography carries it through the stochastic bottleneck of being rare. It is the
    finite-population complement to comparing Lyapunov exponents, which only say
    whether the EXPECTED growth rate is favourable, not whether a lineage of size 1
    survives long enough to enjoy it.
    """
    return _run_replicates(
        rng, n_reps, t_max, sigma_AN, sigma_NA, (0, 0, 1, 0, 0),
        success=lambda ct: ct[U_L] + ct[P_L] >= threshold,
        failure=lambda ct: ct[U_L] + ct[P_L] == 0,
        k_switch=k_switch, q_A_S=q_A_S, q_A_L=q_A_L, q_N=q_N,
        g_A=g_A, g_N=g_N, c=c, D=D, lambda_A=lambda_A, lambda_N=lambda_N,
        delta=0.0, beta=0.0, m=1.0, d0=d0)


def phage_persistence_probability(rng, sigma_AN, sigma_NA, k_switch, q_A_S, q_A_L,
                                   q_N, g_A, g_N, c, D, lambda_A, lambda_N,
                                   delta, beta, m, S0=300, d0=0.5, threshold=40,
                                   n_reps=300, t_max=200.0):
    """P(a single introduced lysogen + its phage persist, reaching `threshold`).

    Runs against a FINITE susceptible background of S0 cells, which -- unlike the
    rare-invasion linearisation in lib.fitness.phage_invasion_growth_rate -- can
    actually be depleted by infection.
    """
    return _run_replicates(
        rng, n_reps, t_max, sigma_AN, sigma_NA, (S0, 0, 1, 0, 0),
        success=lambda ct: ct[U_L] + ct[P_L] + ct[PHAGE] >= threshold,
        failure=lambda ct: ct[U_L] + ct[P_L] + ct[PHAGE] == 0,
        k_switch=k_switch, q_A_S=q_A_S, q_A_L=q_A_L, q_N=q_N,
        g_A=g_A, g_N=g_N, c=c, D=D, lambda_A=lambda_A, lambda_N=lambda_N,
        delta=delta, beta=beta, m=m, d0=d0)


def phage_fitness_decomposition(growth_advantage, lambda_bar, m, delta, S_star, beta):
    """Split the phage's rare-invasion rate into the part it would have with no
    horizontal transmission (its host lineage's relative advantage minus the lysis
    tax) and the residual contribution from infecting new hosts.

    Deliberately a decomposition rather than a weighted "score": choosing weights
    for vertical vs. horizontal transmission would be an invented parameter the
    data does not determine (NOTES.md, open design choice #3).
    """
    vertical_only = growth_advantage - lambda_bar
    total = phage_invasion_growth_rate(growth_advantage, lambda_bar, m, delta,
                                        S_star, beta)
    return {
        "vertical_only": vertical_only,
        "total": total,
        "horizontal_contribution": total - vertical_only,
    }


def _sanity_checks():
    rng = np.random.default_rng(11)
    common = dict(sigma_AN=0.1, sigma_NA=0.1, k_switch=1.0, q_A_S=0.1, q_A_L=0.0,
                  q_N=0.95, g_A=1.0, g_N=0.4, c=0.1, D=5.0)

    # Establishment responds to the lysis tax rather than saturating at 1.0 -- the
    # regression guard for the old pure-birth formulation (see dynamics_stochastic).
    p_free = lysogen_establishment_probability(rng, **common, lambda_A=0.0,
                                                lambda_N=0.0, d0=1.0, n_reps=200)
    p_taxed = lysogen_establishment_probability(rng, **common, lambda_A=0.6,
                                                 lambda_N=0.6, d0=1.0, n_reps=200)
    assert p_taxed < p_free, (p_taxed, p_free)
    assert 0.0 < p_free < 1.0, p_free

    # Horizontal transmission helps a phage that would otherwise struggle.
    kw = dict(**common, lambda_A=0.5, lambda_N=0.5, m=0.5, d0=1.0,
              n_reps=150, t_max=80.0)
    p_with = phage_persistence_probability(rng, delta=0.01, beta=25.0, **kw)
    p_without = phage_persistence_probability(rng, delta=0.0, beta=0.0, **kw)
    assert p_with > p_without, (p_with, p_without)

    # With delta=0 the horizontal contribution is exactly zero.
    d = phage_fitness_decomposition(0.05, 0.02, m=0.5, delta=0.0, S_star=50.0, beta=20.0)
    assert abs(d["horizontal_contribution"]) < 1e-12, d
    assert abs(d["vertical_only"] - 0.03) < 1e-12, d

    print("All sanity checks passed.")


if __name__ == "__main__":
    _sanity_checks()
