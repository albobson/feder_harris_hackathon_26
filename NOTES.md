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
  `k_switch` is an explicit parameter and both genotypes carry it. There is no
  closed form when `k_switch` and `σ` are comparable (a top-Lyapunov-exponent-of-a-
  random-matrix-product problem), so `lib/fitness.py`'s numeric estimator is the
  primary tool, not a validation step. Both genotypes must be evaluated on the
  **same** realized environment path — they compete in one shared environment.

## Both genotypes get the SAME machinery — this is the model's core commitment

Lysogen and non-lysogen are the same object: a population split between unprepared
and prepared cells relaxing at the same `k_switch`. They differ only in the aerobic
target prepared fraction (`q_A` = 0 for the lysogen, ~0.1 for the non-lysogen) and
in the lysogen's lysis rate. Anaerobic targets are identical, because the papers
find anaerobic torCAD expression experimentally indistinguishable between them.

**Do not "simplify" the lysogen back into a single compartment that tracks oxygen
directly.** An earlier draft did exactly that, which meant the lysogen never paid
the post-transition stall penalty the non-lysogen paid — backwards relative to
Carey Fig 1C/D, where it is the *lysogens* that fail to grow after rapid oxygen
depletion. The asymmetry inverted the model's headline result and produced a
"surprising" conclusion that was written into SUMMARY.md before being caught. See
`derivation.md` §4(c).

## Stochastic simulator: no carrying capacity, and turnover matters

`sim/dynamics_stochastic.py` has no population cap. Any compartment with sustained
positive net growth and no depletion will grow without bound, and since it is an
exact Gillespie simulator that means one event per birth — computationally
infeasible long before it is numerically wrong. Always pass an `extinction_check`
covering **both** absorbing outcomes (gone, and past threshold), and keep horizons
short. `sim/dynamics_ode.py` supports an optional `K` on the mean-field side.

Separately: demographic stochasticity depends on birth and death rates
**separately**, not on their difference. The simulator therefore takes a baseline
turnover `d0` and splits a net rate `g` as `b = max(g + d0, 0)`, `d = b − g`. An
earlier version used `max(g,0)` births and `max(−g,0)` deaths — never both — making
every lineage with positive growth a pure-birth process that could not go extinct,
so establishment probability was identically 1.0. `d0` cancels out of growth-rate
comparisons and only affects stochastic quantities.

## Open design choices

1. `D` (stall penalty): currently a growth-rate deduction applied while a cell is
   unprepared under anaerobiosis. The papers describe cells that fail to grow at
   all post-transition, so a discrete stochastic death may be more faithful.
2. Every infection is assumed to produce a lysogen; a real temperate phage
   sometimes goes lytic instead. This makes lysogeny slightly "too easy" for the
   phage.
3. Vertical vs. horizontal transmission are reported as a *decomposition*
   (`sim/metrics.phage_fitness_decomposition`) rather than a weighted score —
   choosing weights would invent a parameter the data does not determine.
4. `q_A` for the lysogen is set to exactly 0. Fig 1B shows a sharp peak at
   background, so this is close, but it is a choice, and the host's entire
   possible benefit lives in the small difference `q_A_S − q_A_L`.
5. Two-state environment with exponential dwell times. `sim/environment.py`
   supports Gamma-distributed dwells (`gamma_shape`) as a separate
   "predictability" axis, which no result currently uses.

## Sensitivity worth knowing before quoting any number

The host's maximum growth advantage from losing the hedge is small — of order
`p_A · q_A · c`, which with the illustrative values (`p_A`≈0.5, `q_A`=0.1,
`c`=0.1) is ~0.005, i.e. half a percent of a generation. Any lysis rate above
about that value erases it entirely. So conclusions about *the host* benefiting
are sensitive to `c` and `q_A` in a way that conclusions about the phage are not.
If the real `c` is much larger than we assumed, the host-benefit window widens; if
it is genuinely undetectable (as the 2018 competition assays suggest), it closes.
