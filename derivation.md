# Derivation: host and phage fitness under environmental switching

This derives the quantities used in `lib/fitness.py`. **This file previously
contained a modeling error, caught during implementation by checking the formulas
against direct simulation — see "A corrected error" below before reading the rest.**
All results below were checked against direct simulation (matrix-product / ODE
integration, not just symbolic algebra) before being finalized.

Notation: environment `E(t) ∈ {A (aerobic), N (anaerobic)}`, a continuous-time Markov
chain with rates `σ_AN` (A→N) and `σ_NA` (N→A), stationary occupancies `p_A`, `p_N`.
Baseline growth rates `g_A`, `g_N` (aerobic growth is generally higher, since
anaerobic TMAO respiration has a lower energy yield). Cost of needless aerobic
torCAD expression: `c`. Penalty for being caught unprepared at the moment
anaerobiosis hits (no fermentable carbon, no torCAD): `D`.

## A corrected error: environment is shared, not a per-cell compartment

The first version of this derivation treated environment the same way it treated
phenotype: as a compartment a population could be distributed across, building a
combined generator over `{(A,U), (A,P), (N,U), (N,P)}` and taking its top eigenvalue.
That's wrong. Oxygen level is **one shared condition** — the whole population is
aerobic together or anaerobic together, never split between the two at a given
moment. Phenotype (prepared/unprepared) genuinely is a per-cell attribute cells can
be distributed across; environment isn't. Mixing the two into one generator silently
assumed each cell had its own independent copy of the environment-switching process,
which isn't the biology.

This was caught by the cheapest possible check: for the lysogen (no phenotype
compartment at all — its torCAD state is a deterministic function of environment,
so population size obeys `dN/dt = r(E(t))·N`), `log(N(t))` is *exactly* the
time-integral of `r(E(s))` — freshman calculus, no eigenvalue possible. Direct
simulation gave a long-run growth rate of 0.648 against a simple weighted average's
prediction of 0.650, while the (wrong) eigenvalue formula gave 0.811. The two
non-eigenvalue numbers agree; the eigenvalue one is just wrong.

## 1. Lysogen: exact result for any parameters

No internal phenotype compartment, no lag: `dN/dt = r(E(t))·N` with
`r(A) = g_A − λ_A`, `r(N) = g_N − λ_N` (`λ` = lysis rate, possibly
environment-dependent — the exploratory oxygen-dependent-lysis hypothesis). By the
ergodic theorem, the long-run growth rate is exactly the stationary time-average:

```
Λ_L = p_A·(g_A − λ_A) + p_N·(g_N − λ_N)
```

Exact for **any** `σ_AN`, `σ_NA` — no timescale assumption, no matrix.

**Verification:** matches a direct path-integral of `r(E(s))` over a long simulated
environment trajectory to within simulation noise (checked in `lib/fitness.py`'s
`_sanity_checks`).

## 2. Non-lysogen: two closed-form limits, numerical in general

The non-lysogen has a genuine per-cell phenotype `X(t) ∈ {U, P}` that relaxes toward
an environment-dependent stationary "prepared fraction" `q_E` at rate `k_switch`
(~1/20 min, per the study's author — not assumed fast or slow relative to `σ`, see
NOTES.md). Within a fixed environment, this is a real 2-compartment linear system
(`v = (U_S, P_S)`, `dv/dt = A(E)·v`) — but `E(t)` is exogenous (shared, not part of
the compartment structure), so the long-run growth rate of `‖v(t)‖` for a randomly
alternating sequence of environment epochs is a genuine **top Lyapunov exponent of a
product of random matrices**. This is a known-hard problem in general (it's exactly
why classical bet-hedging theory, e.g. Kussell & Leibler 2005, restricts its clean
results to limiting regimes) — there is no simple closed form for arbitrary
`k_switch` vs. `σ`.

Two limits **are** closed-form:

- **Quasi-static environment ("fast switching")**: each environment dwell is long
  relative to the (U_S,P_S) system's OWN relaxation time — note this is not simply
  `k_switch ≫ σ`; the relevant relaxation time depends on `k_switch` **and** `c`,
  `D`, `g_A`, `g_N` together (see "a second, subtler bug" below). When it holds,
  the population's growth rate within a dwell settles onto the **top eigenvalue of
  the local (frozen-environment) 2×2 phenotype matrix** `A(state)`, then combines
  across environments by the same plain time-average used elsewhere:
  `Λ_S(fast) = p_A·λtop(A(A)) + p_N·λtop(A(N))`, each `λtop` from the ordinary 2×2
  quadratic formula. This reduces to a plain q-weighted average of the bare
  compartment rates only in the further sub-limit `k_switch → ∞`.
- **Slow switching (`k_switch ≪ σ`)**: phenotype is nearly a fixed lineage trait
  across many environment cycles — the classical bet-hedging regime most of the
  literature assumes (see NOTES.md's fast-switching caveat: this is probably **not**
  the relevant regime here). No local-matrix subtlety here (no mixing happens
  within an epoch in this limit), so it stays a plain weighted average of each
  committed sub-lineage's environment-averaged rate:
  `Λ_S(slow) = max(p_A·g_A + p_N·(g_N−D),  p_A·(g_A−c) + p_N·g_N)`.

**A second, subtler bug, also caught by checking against simulation:** the first
version of the fast-switching formula used a naive q-weighted *arithmetic mean* of
the two compartments' bare growth rates (`g_A − q_A·c`, etc.) instead of the local
matrix's top eigenvalue. That's the correct `k_switch → ∞` asymptote, but a poor
approximation at any finite `k_switch` unless `k_switch` is ALSO large relative to
`c`, `D`, `g_A`, `g_N` — not just large relative to `σ`. Concretely, at
`k_switch=1`, `D=20`, `σ=0.001` (so `k_switch/σ = 1000`, deep in "fast" territory
by the `k_switch≫σ` criterion alone), the naive formula predicted `0.185` while
direct simulation gave `0.665` — off by more than 2×. The corrected (local
eigenvalue) formula gives `0.6645`, matching simulation to within 0.001. This is
the same lesson as the NOTES.md fast-switching caveat, sharpened: "is `k_switch`
fast enough" has to be judged against every rate in the problem, not against `σ`
alone.

For **general `k_switch`** (comparable to `σ` — plausibly exactly the regime that
matters for the free-living/host-lifestyle question), neither limit applies and
there is no closed form. `lib/fitness.py`'s `nonlysogen_growth_rate_numeric`
estimates it directly: propagate a population vector through the exact matrix
exponential of each realized environment segment, renormalizing at every step
(avoids overflow) and accumulating the log-growth. This is mathematically the same
quantity `sim/dynamics_ode.py` computes by brute-force ODE integration (with
`delta=beta=0`) — the closed-form-free `nonlysogen_growth_rate_numeric` is just a
faster, non-stochastic way to get the same number for the pure host-comparison
question, without needing the nonlinear phage machinery.

**Verification:** confirmed the fast-limit formula matches
`nonlysogen_growth_rate_numeric` at `k_switch=10⁴` (to within 0.02) and does *not*
match at `k_switch=10⁻⁴` (confirming the limit is actually limiting, not
accidentally exact everywhere); confirmed the slow-limit formula matches at
`k_switch=10⁻⁴` symmetrically.

## 3. Phage invasion condition (mean-field, linearized, near an all-susceptible population)

Near a population dominated by susceptible non-lysogens at density `S*`, consider a
small lysogen population `L` and free phage `P`. This calculation uses **mean
(time-averaged) rates** throughout — `r_L` = the lysogen's own net growth rate from
torCAD dynamics alone (`lysogen_growth_rate` with `λ_A=λ_N=0`), `λ̄` = the
time-averaged lysis rate (`weighted_time_average(σ_AN, σ_NA, λ_A, λ_N)`) — the same
quasi-static-in-environment spirit as the fast-limit formulas above, not a further
approximation on top of them. Lysogens lyse at rate `λ̄`, releasing burst size `β`;
free phage infect susceptibles at rate `δ·S*` per phage and decay/adsorb at rate `m`.

```
d/dt [L]   [ r_L − λ̄        δ·S*      ] [L]
     [P] = [  β·λ̄      −(m + δ·S*)   ] [P]
```

This is a Metzler matrix (nonnegative off-diagonals): its top eigenvalue is real and
positive **iff `det < 0`**. Expanding and rearranging:

```
Invasion condition:   r_L·(m + δ·S*)  >  λ̄·[m − δ·S*·(β − 1)]
```

Interpretation: the phage-carrying lineage grows when the host's own torCAD-driven
growth advantage (`r_L`) outweighs the lysis tax — but the tax term can go
*negative* when burst size `β` and target density `S*` are large enough
(`δS*(β−1) > m`), meaning horizontal transmission alone can carry a
net-growth-disadvantaged lysogen lineage. This is exactly the host/phage-conflict
structure the model needs to expose: `r_L` (host-favoring) and `β, δ, S*`
(phage-favoring) can pull in different directions.

**Verification:** confirmed the `det<0 ⟺ top eigenvalue>0` equivalence holds over
20,000 random positive parameter draws (0 mismatches); confirmed the `δ=0` (no
horizontal transmission) limit collapses exactly onto `r_L − λ̄`.

## What this does and doesn't settle

- Lysogen fitness (`Λ_L`) is exact for any parameters — a plain weighted average,
  not an approximation.
- Non-lysogen fitness (`Λ_S`) is closed-form only in the fast/slow `k_switch`
  limits; the general (comparable-timescale) case requires the numeric Lyapunov
  estimator or direct simulation — genuinely, not just as an extra validation step.
- Phage invasion condition is a mean-field, rare-phage linearization; it says
  nothing about the interior coexistence state (`S`, `L`, `P` all non-negligible),
  where the mass-action `δ·P·S` term is genuinely nonlinear, and it treats
  environment-averaged rates as constants rather than modeling the (L,P) system's
  own response to environment switching explicitly. That's what the full simulation
  layer (`sim/dynamics_ode.py`, `sim/dynamics_stochastic.py`) is for.
- Structural assumptions (e.g., growth-rate-deduction `D` rather than a discrete
  stochastic death event; two-state environment rather than Gamma-dwell) are open
  design choices, tracked in NOTES.md, not resolved here.
