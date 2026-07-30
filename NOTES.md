# Modeling notes

Running log of assumptions, data gaps, and open questions for the bet-hedging vs.
deterministic-switching (torCAD/HK022) model. See `../[plan]` conversation for full
background; this file tracks what's *established fact* vs. *hypothesis* vs. *free
parameter*, so we don't accidentally treat a swept parameter as a measured one six
months from now.

## Established (from Carey et al. 2018 Cell / 2019 eLife)

- Non-lysogens: aerobic torCAD expression is noisy/bimodal across the population
  (bet-hedging); anaerobic expression is uniformly high (low noise). Oxygen controls
  the *variance* of expression, not really the mean.
- Lysogens (HK022 integrated between torS/torT): a phage-encoded promoter drives torS
  independent of IscR/oxygen, producing excess TorS that keeps torCAD uniformly OFF
  aerobically; anaerobically, torT derepression restores balance and torCAD is
  uniformly ON (same as non-lysogen). Oxygen now controls the *mean*, with low noise
  in either state.
- Post-transition growth on TMAO correlates with pre-transition torCAD expression,
  with an apparent threshold (direct evidence for the benefit-when-needed side).
- ~5% of complete sequenced *E. coli* genomes carry a prophage at this exact
  integration site, across diverse phylogenetic groups (recurrent, convergent).

## Complementary genomic survey (Jeffrey Carey, `jeffrey/SUMMARY.md`)

Jeffrey (paper co-author) ran a comparative-genomics scan across 9 bacterial
species (4 strains each) for the general architecture behind the HK022/torS/torT
case: a normally-short intergenic spacer between a divergently-transcribed gene
pair that's dramatically expanded (>5x, >800bp) by an inserted element in some
strains but not others. Key results, relevant to this model:

- Candidate disruptions turned up **only in *E. coli*** in this sample (12
  candidates across 3 non-reference strains); zero in the other 8 species — but
  the sample is tiny (3 non-reference strains/species), so absence elsewhere is
  not evidence of absence, just not yet observed.
- **Not all disruptions are phage.** Two of the strongest *E. coli* hits
  (*pepE*/*rluF*, *flk*/*pdxB*) are intact prophages, structurally direct analogs
  of HK022 — but two others (*tonB*/*yciI*, *wrbA*/*ymdF*) are IS-element/
  transposon insertions (`istA`/`istB` + recombinases), not phage at all.
- **This matters for the host-vs-phage-benefit question directly**: an IS element
  disrupting the same regulatory architecture would plausibly give the SAME
  host-side benefit (losing costly bet-hedging in a stable niche) WITHOUT the
  lysis risk a prophage carries. If such non-phage routes to the same phenotype
  are similarly accessible mutationally, that's a reason to expect the *phage*
  route to be disfavored for the host specifically wherever an IS-mediated
  alternative is available — which would argue the phage-mediated cases persist
  more because of the *phage's* own interest (creating/maintaining its own
  integration site) than because bacteria are converging on phage lysogeny per se
  as the way to solve this. Not yet built into the model; worth flagging as a
  motivation for eventually adding a third, non-lytic "genotype" (same
  deterministic-switching benefit as the lysogen, but `lambda=0`, no lysis risk)
  to see how it changes which strategy wins.

## Data gaps (confirmed absent from both papers — do not treat as known)

- **Fitness cost of aerobic torCAD expression when not subsequently needed.** The
  2018 Cell paper explicitly tried and failed to detect this in direct competition
  assays. Model parameter `c` — always swept, never fit.
- **Real-world frequency/predictability of aerobic↔anaerobic transitions** in the
  niches these bacteria occupy. Only qualitative speculation exists (mammalian
  urinary tract / animal latrines / marine mammals as TMAO-rich, O2-stable niches).
  Model parameter `σ` (environmental switching rate) — always swept.
- **HK022 induction/lysis rate, and whether it depends on oxygen/environment at
  all.** Not measured for this system. Model parameters `λ0`, `α` — always swept,
  `α` explored over both signs since we have no prior on direction.

## Hypotheses being tested by the model (not established facts)

- Host-driven framing: in-host (low `σ`) niches favor deterministic on/off
  switching because the cost of needless aerobic expression outweighs the rarely-
  needed bet-hedging benefit; free-living (high `σ`) niches favor bet-hedging.
- Phage-driven framing (the paper's own Discussion): HK022 disables host
  bet-hedging to reduce a metabolic burden on the host, which is framed as
  benefiting *phage* replication, not necessarily the host.
- We're treating these as two separate, possibly-conflicting fitness questions
  (host invasion fitness vs. phage R0), not one combined score — it's genuinely
  unclear who benefits, and that's part of what the model should reveal.

## Key modeling constraint: fast phenotype switching

- Per the study's author, torS/torT-driven expression state switching happens on a
  ~20 minute timescale — roughly the bacterial doubling time.
- Classic bet-hedging theory (Kussell & Leibler 2005; most of the ecological
  bet-hedging literature — seed dormancy, persister cells) assumes phenotypes are
  long-lived/heritable across several generations. A ~20 min switching rate is the
  *opposite* regime (comparable to one generation, not several).
- Consequence: we do not assume a fast- or slow-switching quasi-equilibrium.
  `k_switch` is an explicit, swept parameter. See `derivation.md` for the math and
  a real modeling error caught during implementation (an earlier draft treated
  environment as a per-cell compartment cells could be distributed across, like
  phenotype — wrong, since oxygen level is one shared condition, not something
  cells experience independently of each other). The corrected picture: the
  lysogen's growth rate is exact for any `σ` (no internal phenotype compartment,
  so it's a plain weighted time-average). The non-lysogen's growth rate has closed
  forms only in the fast (`k_switch ≫ σ`) and slow (`k_switch ≪ σ`) limits — for
  `k_switch` comparable to `σ` (plausibly the free-living/high-`σ` regime the
  host-lifestyle hypothesis cares about most), there is genuinely no closed form
  (it's a top-Lyapunov-exponent-of-a-random-matrix-product problem), and either
  `lib/fitness.py`'s numeric estimator or full simulation is required — not
  optional, not just a validation step.

## Stochastic simulator limitation: no carrying capacity

`sim/dynamics_stochastic.py` has no population cap. Any compartment left with a
sustained positive net growth rate and no depletion mechanism (e.g. a susceptible
background population when `delta=0`, or any run over a long horizon) will grow
without bound — and since it's an exact Gillespie simulator, that means simulating
one discrete event per birth, which becomes computationally infeasible long before
it becomes numerically wrong. Keep stochastic-simulator horizons short and/or add a
carrying-capacity cap (open design choice, not yet implemented) before running
longer sweeps. `sim/dynamics_ode.py` already supports an optional `K` for exactly
this reason on the mean-field side.

## Open design choices (see plan, "Open design choices" section)

1. `B` (stall penalty): smooth growth-rate deduction vs. true stochastic
   death/extinction event. Papers describe cells that fail to grow at all
   post-transition — leaning toward the latter as more faithful.
2. Whether to track free phage `P(t)` explicitly or fold it into an implicit
   quasi-steady-state infection rate (lighter default; upgrade only if the
   phage-invasion dynamics turn out to be sensitive to it).
3. Weighting of vertical vs. horizontal transmission in the phage fitness score —
   an explicit modeling choice we make, not something the data determines.
