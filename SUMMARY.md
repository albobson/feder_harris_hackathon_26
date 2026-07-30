# Executive Summary: When does bet-hedging lose to an on/off switch?

*Theory-modeling side of the project — companion to Jeffrey's genomic survey in
`jeffrey/SUMMARY.md`. Full math in `derivation.md`; code in `lib/`, `sim/`,
`sweeps/`. No new experimental data went into this. It is a theoretical exploration
meant to generate testable predictions, not a fit to measurements.*

> **Note on a revision.** An earlier version of this summary reported that losing the
> hedge becomes *more* favourable as the environment gets more volatile, and floated
> a fast-cycling experiment to chase it. That was an artifact of a modeling error —
> the lysogen had been given an instant, penalty-free response to oxygen that no cell
> has. With the error fixed the result reverses. The corrected version is below;
> the error is documented in `derivation.md` §4(c). Nothing else in the project
> depended on it.

## The question

Carey et al. (2018, 2019) showed that *E. coli* hedges its bets: under aerobic growth
cells vary randomly in how much TMAO-respiration machinery (torCAD) they carry, so a
minority are always ready if oxygen suddenly vanishes. The HK022 prophage abolishes
this — lysogens are uniformly off aerobically and uniformly on anaerobically.

1. **Would a bacterium ever benefit** from trading the hedge for a clean switch?
2. **Does the phage benefit** — or is this something the phage maintains for its own
   reasons, at the host's expense?

We compute those two things separately and never combine them into a single score,
because they can point in opposite directions, and where they do is the interesting
part.

## What the model is, in plain terms

Both bacteria — with and without the phage — are represented the same way: a
population of cells, each either **prepared** (torCAD on) or **unprepared** (off),
switching between those states on a ~20-minute timescale. Prepared cells pay a small
cost for machinery they aren't using yet. Unprepared cells that get caught by a
sudden loss of oxygen stall until they catch up.

The prophage changes exactly one thing: it drives the aerobic prepared fraction to
**zero**. That's the whole of its modelled effect, and it's exactly what Fig 5
describes. It saves the cost of unnecessary machinery, and it means **nobody is
ready when the oxygen goes** — which is precisely the failure the lysogens show in
Fig 1C/D. Lysogens additionally risk being killed when the prophage activates.

Getting this symmetry right matters. If you let the lysogen respond to oxygen
instantly, you accidentally hand it an advantage no real cell has, and the model
tells you the opposite of the truth. That's the error the note above refers to.

## Finding 1: the window where losing the hedge pays is narrow

![Growth advantage of on/off switching vs. environmental volatility, at three lysis rates](figs/crossover_line.png)

Losing the hedge pays off only when oxygen availability is **stable** — the left-hand
side of the plot, the regime a resident gut or urinary-tract bacterium plausibly
occupies. As flips get more frequent, the advantage falls, crosses zero, and
bet-hedging wins. **This is the direction the original hypothesis predicted**, and it
now falls out of the model rather than being assumed.

Two cautions, both important:

- **The advantage is small even at its best** — under a percent of a generation's
  growth. It is roughly "how often you'd have paid the cost pointlessly," and neither
  that cost nor the payoff has ever been measured in this system.
- **A very small lysis risk erases it.** The green and red curves add lysis rates of
  0.002 and 0.01 per generation. At 0.01 the bacterium never comes out ahead at any
  environmental stability. So whether the *host* benefits hinges on the prophage
  being close to free to carry.

## Finding 2: where the phage spreads, it is usually spreading at the host's expense

![Regions of the volatility / lysis-rate plane showing where the phage spreads and whether the host gains](figs/phase_diagram.png)

Each point is a combination of environmental stability (left–right) and how costly
the prophage is to carry (bottom–top). Three regions:

- **Green, "both gain"** — small, tucked into the stable/very-low-lysis corner. Here
  the bacterium genuinely profits from losing the hedge and the phage spreads too.
- **Red** — the phage spreads while the bacterium is worse off than it would be
  without it. This covers the great majority of the region where the phage can
  spread at all. (Roughly four-fifths of the plotted area, though read that as
  "most of it" rather than a statistic — an area fraction on log axes depends on
  the ranges we chose to plot.)
- **Grey** — the phage cannot establish.

The mechanism behind the red region is worth stating plainly: lysis is how the phage
transmits horizontally. Killing hosts is costly *to hosts* but is the phage's route
to new ones, so a prophage whose host is at a growth disadvantage can still spread —
provided it lyses often enough. The dashed line is exactly that threshold, and it is
why the red region sits *above* it.

This lines up with your comparative-genomics result. Jeffrey found the same
regulatory disruption achieved by plain IS elements and transposons — **no lysis
risk at all**. If a bacterium can get the same rewiring without carrying a phage,
then the phage-mediated cases are hard to explain as something bacteria are selecting
*for*. The model says the same thing from the other direction: the phage does not
need its host to benefit.

## What is assumed vs. measured — please push back

Nothing in the first three below has been measured for this system. We used
illustrative values and checked the qualitative story is not sensitive to them,
except where flagged:

- **Cost of making torCAD when it isn't needed.** The 2018 paper looked for this and
  could not detect it in competition assays, so we used a small value. **The host's
  entire possible benefit is proportional to this number** — it is the single
  parameter that most affects Finding 1.
- **How often oxygen availability actually flips** in a host-associated vs. a
  free-living niche. This is the x-axis of both figures and the crux of the
  in-host/free-living hypothesis. If anyone can bound it even to an order of
  magnitude, that would sharpen this a lot.
- **How often HK022 spontaneously lyses**, and whether that rate depends on oxygen.
  Swept across a wide range; the direction of any oxygen-dependence is unexplored.
- **How fast expression state switches** — this one *is* pinned down, thanks to your
  input: ~20 minutes, about one generation. That is fast enough that standard
  bet-hedging theory's usual assumption (phenotypes stable over many generations)
  does not apply, so we built the model to handle it explicitly rather than assume
  it away.

## Bottom line

Losing the hedge pays for the bacterium only in stable environments, only weakly, and
only if the prophage is almost free to carry — a narrow window. Over most of the
plausible range the phage can spread while leaving its host worse off, because lysis
is simultaneously the host's cost and the phage's means of transmission. Combined
with Jeffrey's finding that IS elements achieve the same rewiring without any lysis
risk, the weight of evidence leans toward **the phage, not the bacterium, being the
beneficiary** of this particular arrangement.

The most useful thing an experiment could add is a bound on the cost of unnecessary
aerobic torCAD expression. It is the parameter the host-side conclusion rests on,
and the one already known to be hard to detect — which is itself informative, since
a cost too small to measure is also a cost too small to make lysogeny worth it.
