# Derivation: host and phage fitness under environmental switching

The math behind `lib/fitness.py`. **Three modeling errors were caught and corrected
during development, each by checking a formula against direct simulation rather than
by re-reading the algebra.** They're documented in §4 because the reasoning that
produced them is easy to repeat.

Notation: environment `E(t) ∈ {A (aerobic), N (anaerobic)}`, a continuous-time Markov
chain with rates `σ_AN`, `σ_NA` and stationary occupancies `p_A`, `p_N`. Growth rates
`g_A`, `g_N` (aerobic growth is higher — TMAO respiration yields less energy). `c` is
the cost of making torCAD before it's needed; `D` the stall penalty for being caught
without it when the oxygen goes; `k_switch` the phenotype relaxation rate (~1/20 min,
about one generation); `λ(E)` the lysis rate.

## 1. One structure, two genotypes

Both genotypes are the same object: a population split between **unprepared** (`U`,
torCAD off) and **prepared** (`P`, torCAD on) cells, relaxing toward an
environment-dependent target prepared fraction `q_E` at rate `k_switch`. Within a
fixed environment,

```
d/dt (U, P) = A(E)·(U, P),   A(E) = [ g_U − λ − k·q      k·(1−q)   ]
                                    [ k·q            g_P − λ − k·(1−q) ]
```

with `(g_U, g_P) = (g_A, g_A − c)` aerobically and `(g_N − D, g_N)` anaerobically.
Growth rates depend on **phenotype and environment only — never on genotype.**

The genotypes differ in exactly two places:

| | non-lysogen | lysogen |
|---|---|---|
| aerobic target `q_A` | ~0.1 (noisy bet-hedging) | **0** (uniformly off) |
| anaerobic target `q_N` | ~0.95 | ~0.95 (same) |
| lysis rate `λ` | 0 | `λ(E)` |

`q_N` is shared because the papers find lysogen and non-lysogen anaerobic torCAD
expression experimentally indistinguishable. So **the prophage's entire modelled
effect is to set the aerobic prepared fraction to zero** — which is exactly what
Carey et al. 2019 Fig 5 describes, and it buys the saved cost `c` at the price of
having nobody ready when oxygen disappears.

## 2. Long-run growth rate

`E(t)` is exogenous and **shared** — one oxygen level for the whole population, not a
compartment cells are distributed across. So the growth of `‖(U,P)‖` along a realized
environment path is a **top Lyapunov exponent of a product of random matrices**: no
closed form in general (precisely why classical bet-hedging theory restricts itself
to limiting regimes). `growth_rate_numeric` estimates it directly by propagating a
population vector through `expm(A(E)·dt)` per segment, renormalizing each step.

Both genotypes must be evaluated on the **same** realized path — they compete in one
shared environment, so the comparison is between their growth along one common
realization, not between separate random draws.

Two closed-form limits, both verified against the numeric estimator:

- **Quasi-static environment**: each dwell is long relative to the phenotype
  system's own relaxation time, so within a dwell growth settles onto
  `λ_top(A(E))`, and `Λ = p_A·λ_top(A(A)) + p_N·λ_top(A(N))`. The validity condition
  is *not* simply `k_switch ≫ σ` — the local relaxation time depends on `k_switch`
  together with `c`, `D`, `g_A`, `g_N`.
- **Frozen phenotype** (`k_switch → 0`): the classical bet-hedging regime most of the
  literature assumes, and per NOTES.md probably *not* the relevant one here. The
  population is dominated by the better committed sub-lineage:
  `Λ = max(p_A·g_A + p_N·(g_N−D), p_A·(g_A−c) + p_N·g_N)`.

Because lysis hits both phenotypes equally it is a uniform tax: `Λ(λ) = Λ(0) − λ̄`
exactly, where `λ̄ = p_A·λ_A + p_N·λ_N`. (Asserted as a test.)

## 3. Phage invasion, measured against the resident

A rare lysogen plus free phage against a resident non-lysogen population at its
ecological equilibrium density `S*`:

```
d/dt [L]   [ s_L − λ̄        δ·S*      ] [L]         s_L = Λ_L(λ=0) − Λ_S
     [P] = [  β·λ̄      −(m + δ·S*)   ] [P]
```

`s_L` is the lysogen's growth advantage **relative to the resident**. In a
density-regulated population the resident's own net per-capita growth is zero at
equilibrium, so a rare lineage spreads iff it beats the resident — not iff it beats
zero. Metzler matrix ⇒ top eigenvalue positive iff `det < 0`, giving

```
s_L·(m + δ·S*)  >  λ̄·[m − δ·S*·(β − 1)]
```

Dividing through by `(m + δS*)` puts it in a more readable form. With
`f = δS*/(m + δS*)` the fraction of released phage that find a host rather than
decaying, and `R₀ = β·f` the successful new lysogens per lysis event:

```
s_L  >  λ̄·(1 − R₀)
```

So when `R₀ > 1` the right side is negative and the phage spreads **even while
harming its host** — and, more subtly, a phage whose host is at a disadvantage
(`s_L < 0`) requires a *minimum* lysis rate `λ̄ > s_L/(1 − R₀)` to invade, because
lysis is how it transmits horizontally. That threshold is the dashed curve in
`figs/phase_diagram.png`.

## 4. The three corrected errors

**(a) Environment treated as a per-cell compartment.** The first draft built a
4-state generator over `{(A,U), (A,P), (N,U), (N,P)}` and took its top eigenvalue.
That treats each cell as having its own private copy of the environment-switching
process, and computes the *annealed* growth rate (the rate of `E[N]` averaged over
environment realizations) rather than the *quenched* one (the a.s. rate of
`log N` along one shared path). For competing two strategies in one environment the
quenched rate is what matters. Caught by noting that a genotype with no phenotype
compartment obeys `dN/dt = r(E(t))·N`, whose log-growth is a plain time-integral —
no eigenvalue is possible. Simulation: 0.648, time-average 0.650, bad formula 0.811.

**(b) Naive average instead of the local eigenvalue in the fast limit.** The
quasi-static formula originally used a `q`-weighted arithmetic mean of the bare
compartment rates, which is only correct as `k_switch → ∞`. At `k_switch = 1`,
`D = 20`, `σ = 0.001` — deep in "fast" territory by the `k_switch ≫ σ` criterion —
it predicted 0.185 against a simulated 0.665. The local-eigenvalue form gives 0.6645.

**(c) The lysogen was given a lag-free, penalty-free response.** This was the
consequential one. The lysogen was modelled as a single compartment tracking oxygen
instantaneously, so it never paid `D`, while the non-lysogen paid `D` on its
unprepared cells throughout every post-transition relaxation. That is backwards
relative to the experiment: Carey Fig 1C/D shows it is the **lysogens** that fail to
grow after rapid oxygen depletion, because nobody is prepared. The asymmetry handed
the lysogen an advantage that grew with environmental volatility — an artifact that
was reported in `SUMMARY.md` as a real and "surprising" finding, and nearly became
the basis of an experimental recommendation.

With both genotypes given the same machinery, the sign reverses:

| σ | `host_gap`, before fix | `host_gap`, after fix |
|---|---|---|
| 0.001 | +0.028 | +0.004 |
| 0.1 | +0.093 | −0.012 |
| 1.0 | +0.309 | −0.049 |
| 10 | +0.416 | −0.048 |

The corrected trend is the intuitive one — stable favours the switch, volatile
favours the hedge — and it saturates near `σ ≈ 3` because once the environment flips
much faster than `k_switch` neither genotype can track it at all.

**(d) Invasion judged against zero rather than the resident.** `phage_invades`
originally asked whether the lysogen+phage lineage grew at all, in a model where the
resident was itself growing exponentially at `Λ_S ≈ 0.64`. That put the invasion
boundary roughly an order of magnitude too high in `λ`. Fixed by the resident-
equilibrium formulation in §3, which also resolves the related inconsistency of
holding `S*` fixed while `S` grew without bound.

## 5. What this settles and what it doesn't

- The host-side comparison (`Λ_L` vs `Λ_S`) is a fair, like-for-like contest, exact
  for any parameters given the model's structure.
- The phage invasion condition is exact in the rare-phage limit against a resident at
  equilibrium. It says nothing about the interior coexistence state where `S`, `L`,
  `P` are all appreciable and `δ·P·S` is genuinely nonlinear — that is what
  `sim/dynamics_ode.py` (with a carrying capacity `K`) and `sim/dynamics_stochastic.py`
  are for.
- Structural choices — `D` as a rate penalty rather than a discrete death, a
  two-state environment rather than Gamma-distributed dwell times, every infection
  producing a lysogen — remain open and are tracked in NOTES.md.
