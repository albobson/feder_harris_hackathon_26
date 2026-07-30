"""Exact stochastic (Gillespie/SSA) simulator for the (U_S, P_S, L, P_phage) system,
driven by an exogenous piecewise-constant environment (sim/environment.py).

This exists for exactly the two things the mean-field ODE (dynamics_ode.py) can't
capture: (1) a rare mutant lineage (e.g. L starting at 1 individual) can go extinct
by chance even with a positive deterministic growth rate -- ODEs are deterministic
and can't show that; (2) phage persistence/extinction is inherently discrete and
bursty (lysis events release an integer burst size at discrete times).

Implementation is plain Gillespie (exact SSA), not tau-leaping -- simpler and exact,
at the cost of scaling with total event count. Fine for the modest population sizes
(tens to a few thousand individuals) relevant to invasion-probability and
persistence-probability studies; not intended for large-population phase-diagram
sweeps (use dynamics_ode.py for those).
"""
import numpy as np


def _rates(U_S, P_S, L, P_phage, state, p):
    q = p["q_A"] if state == "A" else p["q_N"]
    g_S_U = p["g_A"] if state == "A" else p["g_N"] - p["D"]
    g_S_P = (p["g_A"] - p["c"]) if state == "A" else p["g_N"]
    g_L = p["g_A"] if state == "A" else p["g_N"]
    lam = p["lambda_A"] if state == "A" else p["lambda_N"]

    S_total = U_S + P_S

    return {
        "US_birth": max(g_S_U, 0.0) * U_S,
        "US_death": max(-g_S_U, 0.0) * U_S,
        "PS_birth": max(g_S_P, 0.0) * P_S,
        "PS_death": max(-g_S_P, 0.0) * P_S,
        "US_to_PS": p["k_switch"] * q * U_S,
        "PS_to_US": p["k_switch"] * (1 - q) * P_S,
        "L_birth": max(g_L, 0.0) * L,
        "L_death": max(-g_L, 0.0) * L,
        "lysis": lam * L,
        "phage_decay": p["m"] * P_phage,
        "infection": p["delta"] * P_phage * S_total,
    }


def _apply_event(event, counts, rng):
    """Apply all events except "lysis", which the caller handles directly (it
    needs the mean burst size beta, not available here)."""
    U_S, P_S, L, P_phage = counts
    if event == "US_birth":
        U_S += 1
    elif event == "US_death":
        U_S -= 1
    elif event == "PS_birth":
        P_S += 1
    elif event == "PS_death":
        P_S -= 1
    elif event == "US_to_PS":
        U_S -= 1
        P_S += 1
    elif event == "PS_to_US":
        P_S -= 1
        U_S += 1
    elif event == "L_birth":
        L += 1
    elif event == "L_death":
        L -= 1
    elif event == "phage_decay":
        P_phage -= 1
    elif event == "infection":
        S_total = U_S + P_S
        if S_total > 0 and rng.random() < U_S / S_total:
            U_S -= 1
        else:
            P_S -= 1
        L += 1
        P_phage -= 1
    return (U_S, P_S, L, P_phage)


def simulate_stochastic(segments, y0_counts, k_switch, q_A, q_N, g_A, g_N, c, D,
                         lambda_A, lambda_N, delta, beta, m, rng,
                         extinction_check=None, max_events=2_000_000):
    """Run the exact SSA across all environment segments.

    y0_counts: (U_S, P_S, L, P_phage) integer counts.
    extinction_check: optional callable(counts) -> bool; if it returns True the
        simulation stops early (e.g. "L and P_phage both zero" for an invasion
        study once host-plus-phage is unambiguously extinct).
    Returns (t_history, counts_history) as a list of floats and a list of
    (U_S,P_S,L,P_phage) tuples, one entry per event (plus the initial state and
    each segment boundary crossed).
    """
    params = dict(k_switch=k_switch, q_A=q_A, q_N=q_N, g_A=g_A, g_N=g_N, c=c, D=D,
                  lambda_A=lambda_A, lambda_N=lambda_N, delta=delta, m=m)

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
            rates = _rates(*counts, state, params)
            total_rate = sum(rates.values())
            if total_rate <= 0:
                break  # nothing can happen for the rest of this segment
            dt = rng.exponential(1.0 / total_rate)
            if t_local + dt > seg_end:
                break
            t_local += dt
            events, weights = zip(*rates.items())
            weights = np.array(weights)
            event = rng.choice(events, p=weights / weights.sum())
            if event == "lysis":
                U_S, P_S, L, P_phage = counts
                burst = rng.poisson(beta)
                counts = (U_S, P_S, L - 1, P_phage + burst)
            else:
                counts = _apply_event(event, counts, rng)
            counts = tuple(max(0, x) for x in counts)
            n_events += 1
            t_history.append(t_local)
            counts_history.append(counts)
        t_history.append(seg_end)
        counts_history.append(counts)

    return t_history, counts_history


def _sanity_checks():
    import sys
    import os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
    from environment import generate_environment_sequence

    rng = np.random.default_rng(3)
    sigma_AN, sigma_NA = 0.05, 0.05
    t_max = 500.0
    segs = generate_environment_sequence(t_max, sigma_AN, sigma_NA, rng)

    # No lysis, no infection: a decent-sized lysogen population's stochastic mean
    # trajectory should track the deterministic exponential growth roughly.
    # Short horizon deliberately -- with no population cap, exponential growth
    # over a long horizon needs an infeasible number of discrete birth events
    # (Gillespie simulates every single one), so keep t_max small here.
    y0 = (0, 0, 200, 0)
    segs_short = generate_environment_sequence(15.0, sigma_AN, sigma_NA, rng)
    t_hist, c_hist = simulate_stochastic(
        segs_short, y0, k_switch=1.0, q_A=0.1, q_N=0.95, g_A=0.5, g_N=0.2, c=0.3, D=5.0,
        lambda_A=0.0, lambda_N=0.0, delta=0.0, beta=0.0, m=1.0, rng=rng)
    L_final = c_hist[-1][2]
    assert L_final > y0[2], "expected net-positive-growth lysogen population to have grown"

    # Extinction check: a doomed population (huge death rate) should go extinct
    # and extinction_check should halt the simulation early.
    rng2 = np.random.default_rng(4)
    y0_doomed = (0, 0, 5, 0)
    t_hist2, c_hist2 = simulate_stochastic(
        segs, y0_doomed, k_switch=1.0, q_A=0.1, q_N=0.0, g_A=-5.0, g_N=-5.0, c=0.0, D=0.0,
        lambda_A=0.0, lambda_N=0.0, delta=0.0, beta=0.0, m=1.0, rng=rng2,
        extinction_check=lambda counts: counts[2] == 0)
    assert c_hist2[-1][2] == 0, "expected doomed population to go extinct"
    assert t_hist2[-1] < t_max, "extinction_check should have stopped the sim early"

    # Rare-mutant invasion probability sanity: starting from L=1 with a strongly
    # positive growth rate and no phage loss, most replicates should establish
    # (reach a modest threshold) rather than go extinct -- pure luck of the first
    # few divisions is the only source of failure. Stop each replicate as soon as
    # the outcome is decided (extinct, or threshold reached) -- letting an
    # established lineage keep growing unboundedly for the full horizon would
    # require simulating an astronomical number of discrete birth events.
    THRESHOLD = 20
    n_reps, n_established = 200, 0
    for i in range(n_reps):
        rngk = np.random.default_rng(100 + i)
        segsk = generate_environment_sequence(200.0, sigma_AN, sigma_NA, rngk)
        t_h, c_h = simulate_stochastic(
            segsk, (0, 0, 1, 0), k_switch=1.0, q_A=0.1, q_N=0.95, g_A=1.0, g_N=0.8,
            c=0.0, D=0.0, lambda_A=0.0, lambda_N=0.0, delta=0.0, beta=0.0, m=1.0,
            rng=rngk, extinction_check=lambda counts: counts[2] == 0 or counts[2] >= THRESHOLD)
        if c_h[-1][2] >= THRESHOLD:
            n_established += 1
    establishment_rate = n_established / n_reps
    assert establishment_rate > 0.5, f"expected majority establishment, got {establishment_rate}"

    print("All sanity checks passed.")


if __name__ == "__main__":
    _sanity_checks()
