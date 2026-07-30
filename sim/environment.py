"""Two-state (Aerobic/Anaerobic) environment process.

Generates a piecewise-constant environment trajectory as a list of
(t_start, t_end, state) segments, rather than a continuous function -- this lets
sim/dynamics_ode.py integrate exactly segment-by-segment (no discontinuity handling
inside the ODE solver) and lets sim/dynamics_stochastic.py look up the current state
cheaply.

States are the strings "A" (aerobic) and "N" (anaerobic).

Dwell times default to exponential (memoryless switching, rate sigma_AN / sigma_NA).
A Gamma-distributed dwell time with shape `gamma_shape` is also supported, as a
second, independent "predictability" axis: for fixed mean dwell time, shape=1 is
memoryless (exponential), shape>1 makes transitions more regular/predictable
(lower variance in dwell time for the same mean), shape<1 makes them burstier.
"""
import numpy as np


def _sample_dwell_time(rng, rate, gamma_shape):
    mean = 1.0 / rate
    if gamma_shape == 1.0:
        return rng.exponential(mean)
    return rng.gamma(shape=gamma_shape, scale=mean / gamma_shape)


def generate_environment_sequence(t_max, sigma_AN, sigma_NA, rng, start="A", gamma_shape=1.0):
    """Return a list of (t_start, t_end, state) segments covering [0, t_max]."""
    segments = []
    t = 0.0
    state = start
    while t < t_max:
        rate = sigma_AN if state == "A" else sigma_NA
        dwell = _sample_dwell_time(rng, rate, gamma_shape)
        t_end = min(t + dwell, t_max)
        segments.append((t, t_end, state))
        t = t_end
        state = "N" if state == "A" else "A"
    return segments


def state_at(t, segments):
    """Look up the environment state at time t via binary search over segments."""
    starts = [s[0] for s in segments]
    idx = np.searchsorted(starts, t, side="right") - 1
    idx = max(0, min(idx, len(segments) - 1))
    return segments[idx][2]


def stationary_fractions(sigma_AN, sigma_NA):
    """Stationary occupancy (p_A, p_N) of the two-state Markov chain."""
    p_A = sigma_NA / (sigma_AN + sigma_NA)
    p_N = sigma_AN / (sigma_AN + sigma_NA)
    return p_A, p_N


def _sanity_checks():
    rng = np.random.default_rng(0)
    sigma_AN, sigma_NA = 0.05, 0.1
    t_max = 20000.0
    segs = generate_environment_sequence(t_max, sigma_AN, sigma_NA, rng)

    # Segments should be contiguous and cover [0, t_max].
    assert segs[0][0] == 0.0
    assert abs(segs[-1][1] - t_max) < 1e-9
    for (a_start, a_end, _), (b_start, _, _) in zip(segs[:-1], segs[1:]):
        assert abs(a_end - b_start) < 1e-9
    for i in range(len(segs) - 1):
        assert segs[i][2] != segs[i + 1][2]  # states alternate

    # Empirical time-in-A fraction should match stationary p_A.
    time_in_A = sum(e - s for s, e, st in segs if st == "A")
    p_A_empirical = time_in_A / t_max
    p_A_theory, _ = stationary_fractions(sigma_AN, sigma_NA)
    assert abs(p_A_empirical - p_A_theory) < 0.02, (p_A_empirical, p_A_theory)

    # state_at should agree with the segment list at sampled points, including
    # exactly at segment boundaries and at t=0 / t=t_max.
    for s, e, st in segs[::37]:
        assert state_at((s + e) / 2, segs) == st
    assert state_at(0.0, segs) == segs[0][2]
    assert state_at(t_max, segs) == segs[-1][2]

    print("All sanity checks passed.")


if __name__ == "__main__":
    _sanity_checks()
