"""Outcome metrics computed from the simulators in dynamics_ode.py /
dynamics_stochastic.py, plus a couple of thin wrappers around lib/fitness.py.

Deliberately keeps host and phage outcomes as SEPARATE metrics throughout (per the
plan: "never merged into one number") -- a caller assembling a phase diagram should
plot both, not a combined score.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "lib"))
from fitness import phage_invasion_growth_rate  # noqa: E402

sys.path.insert(0, os.path.dirname(__file__))
from environment import generate_environment_sequence  # noqa: E402
from dynamics_stochastic import simulate_stochastic  # noqa: E402


def establishment_probability(rng, sigma_AN, sigma_NA, k_switch, q_A, q_N, g_A, g_N,
                                c, D, lambda_A, lambda_N, delta, beta, m,
                                y0_counts, threshold, n_reps, t_max, track_index=2,
                                extinct_indices=(2,)):
    """Fraction of n_reps stochastic replicates in which the tracked compartment
    (default index 2 = L, the lysogen count) reaches `threshold` before every
    compartment in `extinct_indices` hits zero.

    Each replicate gets an independent environment realization (same sigma_AN,
    sigma_NA) and an independent RNG stream derived from `rng`.
    """
    n_established = 0
    for _ in range(n_reps):
        child_rng = np.random.default_rng(rng.integers(0, 2**63 - 1))
        segs = generate_environment_sequence(t_max, sigma_AN, sigma_NA, child_rng)

        def extinction_check(counts, _idx=extinct_indices, _thr=threshold, _tidx=track_index):
            if counts[_tidx] >= _thr:
                return True
            return all(counts[i] == 0 for i in _idx)

        _, counts_hist = simulate_stochastic(
            segs, y0_counts, k_switch, q_A, q_N, g_A, g_N, c, D,
            lambda_A, lambda_N, delta, beta, m, child_rng,
            extinction_check=extinction_check)
        if counts_hist[-1][track_index] >= threshold:
            n_established += 1
    return n_established / n_reps


def host_invasion_probability(rng, sigma_AN, sigma_NA, k_switch, q_A, q_N, g_A, g_N,
                                c, D, threshold=20, n_reps=200, t_max=500.0):
    """Probability a single rare lysogen mutant (L=1, no phage yet) establishes,
    ignoring lysis entirely -- i.e. purely "does the deterministic-switching host
    strategy survive stochastic loss when rare," the finite-population complement
    to lib.fitness.lysogen_growth_rate/nonlysogen_growth_rate_* (which only say
    whether the EXPECTED rate is positive, not whether a rare lineage of size 1
    survives to benefit from it).
    """
    return establishment_probability(
        rng, sigma_AN, sigma_NA, k_switch, q_A, q_N, g_A, g_N, c, D,
        lambda_A=0.0, lambda_N=0.0, delta=0.0, beta=0.0, m=1.0,
        y0_counts=(0, 0, 1, 0), threshold=threshold, n_reps=n_reps, t_max=t_max,
        track_index=2, extinct_indices=(2,))


def phage_persistence_probability(rng, sigma_AN, sigma_NA, k_switch, q_A, q_N, g_A,
                                    g_N, c, D, lambda_A, lambda_N, delta, beta, m,
                                    S0=200, threshold=50, n_reps=200, t_max=500.0):
    """Probability a single lysogen+phage introduction (L=1, one free phage burst
    not yet released) persists (lysogen-or-phage count reaches `threshold`) rather
    than going extinct (L=0 AND P_phage=0), against a background susceptible
    population of size S0 (held as U_S at the environment-appropriate starting
    split; this is finite, so it CAN be depleted -- unlike the mean-field
    rare-invasion linearization in lib.fitness.phage_invasion_growth_rate, which
    assumes an effectively infinite S* reservoir).
    """
    def extinct_indices_check(counts):
        return counts[2] == 0 and counts[3] == 0

    n_established = 0
    for _ in range(n_reps):
        child_rng = np.random.default_rng(rng.integers(0, 2**63 - 1))
        segs = generate_environment_sequence(t_max, sigma_AN, sigma_NA, child_rng)
        y0 = (S0, 0, 1, 0)

        def extinction_check(counts, _thr=threshold):
            if counts[2] + counts[3] >= _thr:
                return True
            return counts[2] == 0 and counts[3] == 0

        _, counts_hist = simulate_stochastic(
            segs, y0, k_switch, q_A, q_N, g_A, g_N, c, D,
            lambda_A, lambda_N, delta, beta, m, child_rng,
            extinction_check=extinction_check)
        L_final, P_final = counts_hist[-1][2], counts_hist[-1][3]
        if L_final + P_final >= threshold:
            n_established += 1
    return n_established / n_reps


def phage_fitness_decomposition(r_L, lambda_bar, m, delta, S_star, beta):
    """Split the phage's rare-invasion growth rate (lib.fitness.
    phage_invasion_growth_rate) into a vertical-only component (what the lysogen
    lineage would do with zero horizontal transmission, i.e. r_L - lambda_bar) and
    the residual horizontal contribution. Not a weighted "score" -- an explicit
    modeling choice was to keep this as a decomposition rather than invent a
    vertical/horizontal weighting parameter (see NOTES.md open design choice #3).
    """
    vertical_only = r_L - lambda_bar
    total = phage_invasion_growth_rate(r_L, lambda_bar, m, delta, S_star, beta)
    return {
        "vertical_only": vertical_only,
        "total": total,
        "horizontal_contribution": total - vertical_only,
    }


def _sanity_checks():
    rng = np.random.default_rng(11)

    # Strongly favorable host growth, no lysis -> high establishment probability.
    p_est = host_invasion_probability(
        rng, sigma_AN=0.1, sigma_NA=0.1, k_switch=1.0, q_A=0.1, q_N=0.95,
        g_A=1.0, g_N=0.8, c=0.0, D=0.0, threshold=20, n_reps=150, t_max=200.0)
    assert p_est > 0.5, f"expected majority establishment, got {p_est}"

    # Strongly doomed host growth -> near-zero establishment probability.
    p_doomed = host_invasion_probability(
        rng, sigma_AN=0.1, sigma_NA=0.1, k_switch=1.0, q_A=0.1, q_N=0.95,
        g_A=-2.0, g_N=-2.0, c=0.0, D=0.0, threshold=20, n_reps=100, t_max=50.0)
    assert p_doomed < 0.05, f"expected near-zero establishment, got {p_doomed}"

    # Phage persistence: with high burst size and decent infection rate, a
    # favorable phage should persist more often than a beta=0 (no horizontal
    # transmission, pure vertical) control with an otherwise-unfavorable host rate.
    p_phage_favorable = phage_persistence_probability(
        rng, sigma_AN=0.1, sigma_NA=0.1, k_switch=1.0, q_A=0.1, q_N=0.95,
        g_A=0.3, g_N=0.2, c=0.0, D=0.0, lambda_A=0.3, lambda_N=0.3,
        delta=0.02, beta=20.0, m=0.5, S0=100, threshold=30, n_reps=150, t_max=100.0)
    # NOTE: delta=0 means the S background is never depleted by infection, and
    # there's no carrying-capacity cap in this stochastic simulator (see NOTES.md)
    # -- a positive-growth S population left to run for t_max=100 would need to
    # simulate trillions of discrete birth events. Short horizon avoids that; the
    # comparison to p_phage_favorable above only needs L+P_phage's fate, which is
    # decided well before then.
    p_phage_no_horizontal = phage_persistence_probability(
        rng, sigma_AN=0.1, sigma_NA=0.1, k_switch=1.0, q_A=0.1, q_N=0.95,
        g_A=0.3, g_N=0.2, c=0.0, D=0.0, lambda_A=0.3, lambda_N=0.3,
        delta=0.0, beta=0.0, m=0.5, S0=100, threshold=30, n_reps=150, t_max=20.0)
    assert p_phage_favorable > p_phage_no_horizontal, (p_phage_favorable, p_phage_no_horizontal)

    # Decomposition: with delta=0 (no horizontal transmission), the horizontal
    # contribution must be exactly zero.
    decomp = phage_fitness_decomposition(r_L=0.2, lambda_bar=0.05, m=0.5, delta=0.0,
                                          S_star=2.0, beta=20.0)
    assert abs(decomp["horizontal_contribution"]) < 1e-10, decomp

    print("All sanity checks passed.")


if __name__ == "__main__":
    _sanity_checks()
