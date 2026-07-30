# Slide deck — modeling side

Four slides plus one backup. Written to be pasted into Google Slides / PowerPoint,
or rendered directly with Marp (`marp SLIDES.md`) or reveal.js — `---` separates
slides. Each slide has a message-style title, a few bullets, and speaker notes.

Audience assumption: experimentalists, no modeling background. Notation is kept out
of the slide bodies on purpose; the notes carry anything technical that might come
up in questions. Full detail lives in `SUMMARY.md` and `derivation.md`.

---

## The model: one knob, two consequences

Every cell is in one of two states, flipping back and forth on a **~20 min**
timescale (about one generation — thanks Jeff):

|  | **aerobic** | **anaerobic** |
|---|---|---|
| **not ready** (torCAD off) | grows fastest | **stalls** — can't respire TMAO yet |
| **ready** (torCAD on) | pays a small cost | grows |

- **Wild-type:** a minority of cells (~10%) stay ready aerobically. That's the hedge.
- **Lysogen:** the prophage drives that fraction to **zero**. That is the *only*
  thing it changes.
- So the phage buys one thing and sells another: it **saves the cost** of unneeded
  machinery, and it leaves **nobody ready** when the oxygen goes.
- We simulate both strains in the *same* fluctuating oxygen environment and ask who
  grows faster over the long run — and, separately, whether the phage itself spreads.

> **Notes.** Both strains get identical machinery and identical switching speed —
> the only differences are the aerobic ready-fraction (0 vs ~10%) and that lysogens
> can be killed when the prophage activates. Anaerobic behaviour is identical in the
> model because your data show it's experimentally indistinguishable.
> This symmetry is the whole ballgame: an earlier version let the lysogen respond to
> oxygen instantly, which quietly gave it a superpower no cell has and produced the
> opposite answer. Reproducing your Fig 1C/D — lysogens failing after rapid O₂
> depletion — is the check that the structure is right.
> If asked "what does 'long-run growth rate' mean": average growth over a long
> stretch of randomly alternating aerobic/anaerobic periods, both strains
> experiencing the *same* sequence.

---

## Losing the hedge only pays in a stable environment — and barely

![](figs/crossover_line.png)

- **Left→right:** how often oxygen flips. Left = stable (a resident gut/urinary-tract
  bacterium). Right = volatile (free-living).
- **Above the grey line** the phage's on/off switch wins; **below it**, bet-hedging wins.
- The crossover **is there, and in the predicted direction** — stable environments
  favour the switch, volatile ones favour the hedge. That came out of the model
  rather than being put in.
- But: the advantage is **under 1% of a generation** even at its best, and the two
  lower curves show that a lysis rate of just **0.01 per generation wipes it out
  entirely** at any environmental stability.

> **Notes.** The three curves are the same calculation at three prophage costs
> (lysis rates 0, 0.002, 0.01). Lysis shifts the curve straight down — it's a flat
> tax, so it can't change the shape, only whether the curve clears zero.
> The dip-then-flatten on the right is real: once oxygen flips much faster than cells
> can switch, neither strain tracks it and the gap saturates.
> Honest caveat to volunteer: the *size* of the advantage is proportional to the cost
> of unneeded torCAD, which nobody has measured — the 2018 competition assays
> couldn't detect it. So treat the y-axis as illustrative and the *sign change* as
> the result.

---

## Where the phage spreads, it's usually spreading at your expense

![](figs/phase_diagram.png)

- Same left→right axis. **Bottom→top:** how costly the prophage is to carry.
- **Green — both gain.** Real, but a small corner: stable environment *and* an almost
  free prophage.
- **Red — the phage spreads while the bacterium is worse off.** The great majority of
  where the phage can spread at all.
- **Grey — the phage can't establish.**
- **Why:** lysis is simultaneously the host's cost *and* the phage's only route to new
  hosts. A prophage whose host is at a disadvantage can still spread — as long as it
  lyses often enough. That's the dashed line, and it's why red sits *above* it.

> **Notes.** The dashed curve is a threshold on how many new lysogens each lysis
> event produces; above ~1 the phage spreads regardless of what it does to the host.
> "Spreads" here means *increases in frequency* — i.e. out-competes the resident
> population, not merely grows.
> Don't over-read the area fractions: it's a log-log plot over ranges we chose.
> The claim is directional — the green region is small and requires special
> conditions — not "78.0% of parameter space."

---

## Putting it together with the genomics: the phage doesn't need you to benefit

**The model says:** the bacterium's benefit window is narrow (stable niche, tiny
cost, near-free prophage). The phage's is wide, and doesn't require host benefit.

**Jeff's survey says:** the same architecture — a mobile element landing in a
divergent gene pair's shared promoter — recurs across *E. coli*, via **more than one
class of element**: two intact prophages (*pepE/rluF*, *flk/pdxB*), plus IS-family
transposase/recombinase insertions (*tonB/yciI*, *wrbA/ymdF*) that are **not intact
prophages** and can't lyse.

**Together:**

- If the *bacterium* were the beneficiary, the favoured route is obviously the one
  **without lysis risk**. Lysis-incompetent elements should dominate these sites.
- Instead, intact prophages persist there — and Carey 2019 found the HK022 site
  occupied in ~5% of sequenced genomes, at a **conserved position across diverse
  phylogroups**, with sequence conservation around the phage promoter.
- That pattern fits **selection on the phage side**, maintaining its own integration
  site, better than it fits bacteria selecting for a phage to do this job.

**Two things that would test it, using data you already have:**

1. At these divergent spacers, is the ratio of **intact prophage : lysis-incompetent
   element** what you'd expect if hosts were driving it? Jeff's 5 unclassified small
   candidates are exactly the cryptic/degraded elements this question turns on.
2. Do prophage-occupied strains skew toward **stable-oxygen niches** (gut, urinary
   tract) relative to IS-occupied ones? Isolation-source metadata is already in the
   Carey 2019 supplement.

**Most valuable single measurement:** a bound on the cost of unneeded aerobic torCAD.
Everything on the host side scales with it — and the fact that it's already known to
be hard to detect is itself informative: a cost too small to measure is a cost too
small to make lysogeny worth it.

> **Notes.** State clearly that the joint inference is an *inference*, not a
> demonstration — the model shows host benefit isn't required, and the genomics shows
> lysis-free alternatives exist; neither alone proves phage-side selection.
> Fair pushback to expect: "your cost parameter is made up." Correct — and the
> response is that the *conclusion* (narrow host window) gets stronger, not weaker,
> as the cost gets smaller, which is the direction the failed competition assays
> point.
> Careful wording: the two IS cases do carry some phage-derived genes (baseplate,
> tape-measure fragments), so they're composite/degraded rather than phage-free.
> The argument only needs them to be lysis-incompetent.

---

## BACKUP — what the stochastic simulations are for

*Not in the main flow: no result in this deck depends on them. Keep in reserve for
"did you account for randomness / small populations?"*

- The main results use **average growth rates** — appropriate when populations are
  large.
- A separate exact **cell-by-cell simulator** covers what averages can't: a lineage
  that starts as **one cell** can be lost by pure chance even when its average growth
  rate is favourable.
- It tracks individual birth, death, switching, infection and lysis events, so it can
  answer: *does a single new lysogen actually establish?* and *does an introduced
  phage persist or die out?*
- Cross-checked against the analytical model, and against textbook birth–death theory
  for establishment probability.

> **Notes.** This layer caught a genuine bug: an earlier version couldn't produce
> extinction at all (it modelled births without deaths), so every lineage "always
> established." Now fixed and validated.
> If pressed on why it's not in the main results: the questions this deck answers are
> about long-run competition in large populations, where averages are the right tool.
> The stochastic layer matters for the *origin* question — whether a single new
> lysogen survives its first few divisions — which we haven't pushed on yet.
