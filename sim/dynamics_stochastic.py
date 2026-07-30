"""Exact stochastic (Gillespie/SSA) simulator for the two-genotype system, driven by
an exogenous piecewise-constant environment (sim/environment.py).

State: (U_S, P_S, U_L, P_L, P_phage) -- unprepared/prepared non-lysogens (S),
unprepared/prepared lysogens (L), and free phage. Both genotypes carry a phenotype
split and relax at the same k_switch; they differ only in their aerobic target
prepared fraction (q_A_S vs q_A_L) and in that lysogens lyse. See lib/fitness.py.

This exists for what the deterministic layer cannot do: a rare lineage starting from
a single cell can go extinct by chance even with a positive growth rate, and phage
persistence is inherently discrete and bursty.

BIRTH/DEATH, and why there is a baseline turnover parameter `d0`
----------------------------------------------------------------
Demographic stochasticity depends on the birth and death rates SEPARATELY, not on
their difference. A pure-birth process (death rate 0) started from one cell can
never go extinct, so its establishment probability is identically 1. An earlier
version of this file derived event rates as `max(g,0)` births and `max(-g,0)` deaths
-- never both at once -- and therefore reported an establishment probability of
exactly 1.0 for every positive growth rate, which is not a result, it is the
absence of one.

Here a net rate g is split as

    b = max(g + d0, 0)        d = b - g

so b - d == g exactly (deterministic behaviour is untouched) while turnover is real.
d0 is a baseline per-capita death rate: raise it for a high-turnover population with
more demographic noise, lower it toward 0 to recover the old degenerate behaviour.
Because it shifts every genotype's net rate identically, d0 cancels out of any
growth-rate COMPARISON -- it only affects stochastic quantities.
"""
import numpy as np

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "lib"))
from fitness import phenotype_growth_rates  # noqa: E402

# State vector indices.
U_S, P_S, U_L, P_L, PHAGE = 0, 1, 2, 3, 4


def _birth_death(g, d0):
    b = max(g + d0, 0.0)
    return b, b - g


def _event_rates(counts, state, p):
    """Returns (names, rates) for every possible event in the current state."""
    g_U, g_P = phenotype_growth_rates(state, p["g_A"], p["g_N"], p["c"], p["D"])
    q_S = p["q_A_S"] if state == "A" else p["q_N"]
    q_L = p["q_A_L"] if state == "A" else p["q_N"]
    lam = p["lambda_A"] if state == "A" else p["lambda_N"]
    d0 = p["d0"]

    bU, dU = _birth_death(g_U, d0)
    bP, dP = _birth_death(g_P, d0)

    n_S = counts[U_S] + counts[P_S]

    names = (
        "US_birth", "US_death", "PS_birth", "PS_death",
        "US_to_PS", "PS_to_US",
        "UL_birth", "UL_death", "PL_birth", "PL_death",
        "UL_to_PL", "PL_to_UL",
        "UL_lysis", "PL_lysis",
        "phage_decay", "infect_US", "infect_PS",
    )
    rates = np.array([
        bU * counts[U_S], dU * counts[U_S],
        bP * counts[P_S], dP * counts[P_S],
        p["k_switch"] * q_S * counts[U_S], p["k_switch"] * (1 - q_S) * counts[P_S],
        bU * counts[U_L], dU * counts[U_L],
        bP * counts[P_L], dP * counts[P_L],
        p["k_switch"] * q_L * counts[U_L], p["k_switch"] * (1 - q_L) * counts[P_L],
        lam * counts[U_L], lam * counts[P_L],
        p["m"] * counts[PHAGE],
        p["delta"] * counts[PHAGE] * counts[U_S],
        p["delta"] * counts[PHAGE] * counts[P_S],
    ], dtype=float)
    return names, rates


def _apply(event, counts, p, rng):
    c = list(counts)
    if event == "US_birth":
        c[U_S] += 1
    elif event == "US_death":
        c[U_S] -= 1
    elif event == "PS_birth":
        c[P_S] += 1
    elif event == "PS_death":
        c[P_S] -= 1
    elif event == "US_to_PS":
        c[U_S] -= 1; c[P_S] += 1
    elif event == "PS_to_US":
        c[P_S] -= 1; c[U_S] += 1
    elif event == "UL_birth":
        c[U_L] += 1
    elif event == "UL_death":
        c[U_L] -= 1
    elif event == "PL_birth":
        c[P_L] += 1
    elif event == "PL_death":
        c[P_L] -= 1
    elif event == "UL_to_PL":
        c[U_L] -= 1; c[P_L] += 1
    elif event == "PL_to_UL":
        c[P_L] -= 1; c[U_L] += 1
    elif event == "UL_lysis":
        c[U_L] -= 1; c[PHAGE] += rng.poisson(p["beta"])
    elif event == "PL_lysis":
        c[P_L] -= 1; c[PHAGE] += rng.poisson(p["beta"])
    elif event == "phage_decay":
        c[PHAGE] -= 1
    elif event == "infect_US":
        # An infected cell keeps its phenotype and becomes a lysogen.
        c[U_S] -= 1; c[U_L] += 1; c[PHAGE] -= 1
    elif event == "infect_PS":
        c[P_S] -= 1; c[P_L] += 1; c[PHAGE] -= 1
    return tuple(max(0, x) for x in c)


def simulate_stochastic(segments, y0_counts, k_switch, q_A_S, q_A_L, q_N,
                         g_A, g_N, c, D, lambda_A, lambda_N, delta, beta, m, rng,
                         d0=0.5, extinction_check=None, max_events=2_000_000):
    """Run the exact SSA across all environment segments.

    y0_counts: (U_S, P_S, U_L, P_L, P_phage) integer counts.
    extinction_check: optional callable(counts) -> bool; stop early when it is True
        (use it for BOTH absorbing outcomes, e.g. "lineage gone OR past threshold",
        otherwise an established lineage grows without bound and the SSA has to
        simulate every single birth -- see NOTES.md on the missing carrying cap).
    Returns (t_history, counts_history).
    """
    p = dict(k_switch=k_switch, q_A_S=q_A_S, q_A_L=q_A_L, q_N=q_N,
             g_A=g_A, g_N=g_N, c=c, D=D, lambda_A=lambda_A, lambda_N=lambda_N,
             delta=delta, beta=beta, m=m, d0=d0)

    counts = tuple(int(x) for x in y0_counts)
    t_history = [segments[0][0]]
    counts_history = [counts]
    n_events = 0

    for seg_start, seg_end, state in segments:
        t_local = seg_start
        while t_local < seg_end:
            if extinction_check is not None and extinction_check(counts):
                return t_history, counts_history
            if n_events >= max_events:
                return t_history, counts_history
            names, rates = _event_rates(counts, state, p)
            total_rate = rates.sum()
            if total_rate <= 0:
                break
            dt = rng.exponential(1.0 / total_rate)
            if t_local + dt > seg_end:
                break
            t_local += dt
            idx = np.searchsorted(np.cumsum(rates), rng.random() * total_rate)
            counts = _apply(names[min(idx, len(names) - 1)], counts, p, rng)
            n_events += 1
            t_history.append(t_local)
            counts_history.append(counts)
        t_history.append(seg_end)
        counts_history.append(counts)

    return t_history, counts_history


def _sanity_checks():
    sys.path.insert(0, os.path.dirname(__file__))
    from environment import generate_environment_sequence

    sigma = 0.1
    g_A, g_N, c, D = 1.0, 0.4, 0.1, 5.0

    def establishment(net_g, d0, n_reps=400, threshold=25, t_max=120.0, seed=0):
        """P(a single lysogen cell founds a lineage reaching `threshold`)."""
        n_ok = 0
        for i in range(n_reps):
            r = np.random.default_rng(seed + i)
            segs = generate_environment_sequence(t_max, sigma, sigma, r)
            _, hist = simulate_stochastic(
                segs, (0, 0, 1, 0, 0), k_switch=1.0, q_A_S=0.1, q_A_L=0.0, q_N=0.95,
                g_A=net_g, g_N=net_g, c=0.0, D=0.0, lambda_A=0.0, lambda_N=0.0,
                delta=0.0, beta=0.0, m=1.0, rng=r, d0=d0,
                extinction_check=lambda ct: (ct[U_L] + ct[P_L]) == 0
                or (ct[U_L] + ct[P_L]) >= threshold)
            n_ok += (hist[-1][U_L] + hist[-1][P_L]) >= threshold
        return n_ok / n_reps

    # THE regression test for the birth/death fix: for a birth-death process the
    # establishment probability from one cell is ~ (b-d)/b, so with real turnover
    # a modestly-advantageous lineage must FAIL a good fraction of the time. The
    # old max(g,0)/max(-g,0) formulation made this identically 1.0.
    p_lo_turnover = establishment(net_g=0.5, d0=0.0)
    p_hi_turnover = establishment(net_g=0.5, d0=2.0)
    assert p_lo_turnover > 0.99, p_lo_turnover          # no deaths -> cannot go extinct
    assert 0.05 < p_hi_turnover < 0.75, p_hi_turnover   # real turnover -> real risk
    expected = 0.5 / (0.5 + 2.0)                        # (b-d)/b for b=2.5, d=2.0
    assert abs(p_hi_turnover - expected) < 0.12, (p_hi_turnover, expected)

    # Establishment must increase with the growth advantage.
    assert establishment(net_g=0.2, d0=2.0) < establishment(net_g=1.5, d0=2.0)

    # Infection moves cells S -> L while conserving total cell number. All growth
    # switched off (g=c=D=d0=0) so the only events are switching and infection.
    r = np.random.default_rng(11)
    segs = generate_environment_sequence(6.0, sigma, sigma, r)
    _, hist = simulate_stochastic(
        segs, (200, 0, 0, 0, 400), k_switch=1.0, q_A_S=0.1, q_A_L=0.0, q_N=0.95,
        g_A=0.0, g_N=0.0, c=0.0, D=0.0, lambda_A=0.0, lambda_N=0.0,
        delta=0.02, beta=0.0, m=0.0, rng=r, d0=0.0)
    final = hist[-1]
    assert final[U_L] + final[P_L] > 0, "expected some infections"
    assert final[U_S] + final[P_S] + final[U_L] + final[P_L] == 200

    # Lysis converts lysogens into free phage.
    r = np.random.default_rng(12)
    segs = generate_environment_sequence(4.0, sigma, sigma, r)
    _, hist = simulate_stochastic(
        segs, (0, 0, 100, 0, 0), k_switch=0.0, q_A_S=0.1, q_A_L=0.0, q_N=0.95,
        g_A=0.0, g_N=0.0, c=0.0, D=0.0, lambda_A=0.5, lambda_N=0.5,
        delta=0.0, beta=10.0, m=0.0, rng=r, d0=0.0)
    assert hist[-1][U_L] < 100 and hist[-1][PHAGE] > 0

    print("All sanity checks passed.")


if __name__ == "__main__":
    _sanity_checks()
