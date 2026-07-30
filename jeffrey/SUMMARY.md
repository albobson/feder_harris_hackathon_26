# Divergent-Promoter Disruption Survey

## Part 1: Cross-strain intergenic-spacer scan across 9 bacterial species

### Motivation

Carey et al. 2019 (eLife) showed that bacteriophage HK022 integrates precisely
within the ~130bp intergenic region shared by the divergently-transcribed
*E. coli* genes *torS* and *torT*. That region normally carries a shared IscR
repressor binding site regulating both genes; HK022 disrupts it and supplies
its own promoter for *torS*, rewiring the cell's TMAO-dependent respiratory
bet-hedging strategy. This raised a broader question: how common is this
mechanism — an inserted mobile element (prophage, transposon, IS element, etc.)
landing in and disrupting the shared regulatory region of a divergent gene
pair — across other bacterial species and other classes of insertion element?

### Method

Insertion elements preferentially target short, AT-rich, protein-free
intergenic stretches — the same property that lets a stretch serve as a
compact shared promoter/operator for a divergent gene pair. That gives a
concrete computational signature to search for: **a normally-short intergenic
spacer between a divergently-transcribed gene pair that is dramatically
expanded, with extra genes crammed into it, in some strains of a species but
not others.**

Pipeline (`analysis/fetch_genomes.py`, `analysis/find_divergent_disruptions.py`):

1. For each species, download several complete RefSeq genome assemblies
   (GFF3 + FASTA) from NCBI's Datasets REST API.
2. Pick one well-known "baseline" strain per species (e.g. *E. coli* K-12
   MG1655, *P. aeruginosa* PAO1, *M. tuberculosis* H37Rv) and, in its
   annotation, find all pairs of **immediately-adjacent genes transcribed
   divergently** (left gene on `-` strand, right gene on `+` strand) with an
   intergenic spacer ≤500bp — candidate shared regulatory regions, matched
   across strains by gene symbol (not locus_tag, which differs per assembly).
3. For every other strain of that species, look up the same gene-symbol pair
   and measure the distance between them. If the pair is no longer adjacent
   and the gap has grown >5x the baseline (and >800bp absolute), flag it as a
   candidate disruption, and record the genes that now sit between them.
4. Scan those intervening genes' product annotations for mobile-element
   keywords (integrase, transposase, phage structural genes, recombinase,
   IS-family names, etc.) to characterize what kind of element is responsible.

**Bug fixed along the way:** distances were initially computed as linear
coordinate differences, which produces spurious multi-megabase "expansions"
when a gene pair happens to sit near wherever a given assembly's arbitrary
circular-chromosome coordinate origin was drawn. Fixed by computing the
shorter of the two directions around the circular replicon (using each
replicon's length from the GFF's `##sequence-region` line). A ceiling of
200kb was also applied to exclude genuine large-scale genome rearrangements
(a >1.6Mb rearrangement in *B. fragilis* and ~3Mb replichore inversions in
*B. subtilis* were caught and correctly excluded this way — real structural
biology, but not a discrete insertion element).

### Data

36 complete RefSeq genomes (4 strains each) spanning 9 species chosen for
breadth across bacterial classes and independent reasons to expect
insertion-element activity:

| Species | Class/phylum | Note |
|---|---|---|
| *Escherichia coli* | Gammaproteobacteria | prophage + IS rich |
| *Pseudomonas aeruginosa* | Gammaproteobacteria | genomic islands |
| *Vibrio cholerae* | Gammaproteobacteria | CTXφ integration |
| *Staphylococcus aureus* | Bacillota (Firmicutes) | φ13/SCCmec |
| *Bacillus subtilis* | Bacillota (Firmicutes) | SPβ, ICEBs1 |
| *Streptococcus pyogenes* | Bacillota (Firmicutes) | prophage virulence factors |
| *Mycobacterium tuberculosis* | Actinobacteria | IS6110-rich genome |
| *Bacteroides fragilis* | Bacteroidota | — |
| *Campylobacter jejuni* | Campylobacterota | — |

Strain selection within each species was a convenience sample: the first 4
distinct-strain "Complete Genome" RefSeq entries returned by the NCBI API,
not a curated or random sample. Data lives in `data/genomes/<species>/<accession>/`;
the manifest of exactly which accessions were used is in `data/genome_manifest.tsv`.

### Results

Candidate disruptions surviving the size filter and orientation checks were
found **only in *E. coli*** (12 candidates, across 3 of the 4 non-reference
strains); the other 8 species showed none in this small sample. Strongest
*E. coli* hits:

| Gene pair | Strain | Baseline spacer | Expanded to | Intervening annotation |
|---|---|---|---|---|
| *pepE*/*rluF* | O157:H7 Sakai (GCF_000008865.2) | 211bp | 45.2kb (62 genes) | full prophage: portal, capsid, tail, tape-measure, baseplate |
| *flk*/*pdxB* | 14EC020 (GCF_002853715.1) | 98bp | 36.3kb (48 genes) | full prophage: integrase, terminase, portal, capsid, tail |
| *tonB*/*yciI* | 14EC020 | 223bp | 14.1kb (14 genes) | IS21-family transposase (*istA*/*istB*) + recombinases + phage baseplate |
| *wrbA*/*ymdF* | 14EC020 | 372bp | 11.4kb (12 genes) | transposase + recombinase + phage tape-measure/baseplate |

Plus 5 smaller *E. coli* candidates (2–7kb insertions) whose intervening genes
are annotated only as hypothetical proteins — not yet classified as
phage/IS/other, but structurally consistent with a cryptic or degraded
mobile element.

The *pepE*/*rluF* and *flk*/*pdxB* cases are structurally direct analogs of
the HK022/*torS*/*torT* mechanism (intact prophage landing in a divergent
gene pair's spacer), while *tonB*/*yciI* and *wrbA*/*ymdF* show the same
outcome from a different class of element (IS/composite transposon rather
than an intact phage) — i.e. the same regulatory-disruption phenomenon is not
specific to phage integration.

### Caveats / limitations

- **Sample size per species is tiny (3 non-reference strains).** A zero-hit
  species does *not* mean the phenomenon is absent there — it means none of
  the 3 sampled strains happened to carry an insertion at that species'
  divergent loci. Species with well-documented high insertion-element
  turnover (*M. tuberculosis*/IS6110, *S. aureus*/φ13, *Shigella*-type IS
  expansion) are the most promising places to add more strains first.
- **Gene-symbol matching** only finds pairs where both flanking genes have an
  assigned symbol in RefSeq annotation; pairs of poorly-characterized genes
  are invisible to this method. A protein-BLAST-based ortholog map (feasible
  with the `blastp`/`makeblastdb` already available locally) would catch more.
- Reference-strain divergent-pair discovery only checks linearly-adjacent
  genes and does not itself check the wraparound pair at each replicon's
  coordinate seam (a minor, at-most-one-pair-per-replicon gap).
- The 5 unclassified small *E. coli* candidates warrant a closer look
  (e.g. against ISfinder/PHASTER) before concluding what they are.

## Part 2: Are the 12 *E. coli* candidate pairs functionally coupled (like *torS*/*torT*)?

### Motivation

In Carey et al. 2019, *torS* and *torT* aren't just genomic neighbors — they
work together in the *same* two-component signaling pathway (TorS is the
sensor kinase, TorT is the periplasmic partner that helps it detect TMAO).
That functional coupling is part of what makes the phage-disruption story
biologically meaningful: the phage isn't just landing between two random
genes, it's rewiring a specific signal-transduction circuit. The natural
follow-up question for the 12 candidate disruptions found in Part 1: do any
of those divergent pairs show the same kind of functional relationship?

### Method

For each of the 9 unique gene pairs underlying the 12 *E. coli* candidates,
a literature/database check (EcoCyc, UniProt, NCBI Gene, primary literature)
was run to look for any documented functional coupling — shared pathway,
shared regulator/regulon, physical complex membership, substrate-product
relationship, or other established genetic/physical interaction — versus no
known connection (i.e., the pair are presumably coincidental genomic
neighbors, and the insertion hotspot reflects spacer/sequence properties
rather than shared regulatory logic).

### Results

| Gene pair | Verdict | Basis |
|---|---|---|
| ybjP / ybjQ | No known link | TolC-lipoprotein export factor vs. uncharacterized UPF0145 protein |
| dusC / yohJ | No known link | tRNA dihydrouridine synthase vs. uncharacterized membrane protein |
| fabR / sthA | No known link | FabR's confirmed regulon (fabA/fabB/yqfA) doesn't include sthA |
| pepE / rluF | No known link | Dipeptidase vs. rRNA/tRNA pseudouridine synthase; no shared regulator/complex |
| dnaB / qorA | No known link | Replicative helicase vs. quinone detoxification enzyme |
| wrbA / ymdF | Weak, unconfirmed | Both are RpoS (general stress regulon) targets, but no *E. coli*-specific study directly links the two genes |
| tonB / yciI | No known link found | TonB energy-transduction complex vs. yciI (function essentially uncharacterized — weaker negative given how poorly annotated yciI is) |
| flk / pdxB | No known link — and independently documented as coincidental | Müller et al. 1997 (PubMed 9079927) explicitly describe this pair sharing a divergent promoter despite completely unrelated functions (flagellar-assembly checkpoint vs. vitamin B6 biosynthesis) |
| btsT / tsr | No known link | Pyruvate:H+ symporter vs. serine chemoreceptor; both loosely relate to cell energy status but no direct regulatory/physical connection found |

### Interpretation

None of the 12 *E. coli* candidates show a functional coupling analogous to
*torS*/*torT*. Where a relationship exists at all (*wrbA*/*ymdF*), it's a
shared stress regulon at best, not a shared pathway or complex. The
*flk*/*pdxB* pair is a useful negative control: it's independently documented
in the literature as two functionally unrelated genes that happen to share a
divergent promoter architecture purely by genomic coincidence — the opposite
of the *torS*/*torT* case.

This suggests *torS*/*torT*'s functional coupling may be the exception rather
than the rule for divergent pairs whose shared regulatory region gets
disrupted by an inserted element: insertion hotspots more plausibly track
spacer/sequence properties (short, AT-rich, protein-free intergenic DNA) than
any shared regulatory logic between the two flanking genes. Finding a true
second example of the *torS*/*torT* phenomenon — a functionally coupled
divergent pair disrupted by an insertion — would likely require scanning many
more strains/species than the small sample used in Part 1.

## Part 3: Expanding to many more genomes per species

### Motivation

Part 1 used only 3 non-reference strains per species — far too few to
conclude a species lacks the phenomenon just because none of the 3 sampled
strains happened to carry an insertion at a given locus. The obvious fix:
sample many more strains per species.

### Method changes

- `fetch_genomes.py` was changed from picking the first 4 "Complete Genome"
  RefSeq records the NCBI API returned (which is not documented as random,
  and risked clustering genomes submitted together by the same study/center)
  to pulling a pool of up to 1000 complete-genome records per species,
  deduplicating by strain, and taking a **random sample of 20 strains** per
  species (fixed random seed for reproducibility). Total: **180 genomes**
  (20 x 9 species) versus 36 before.
- Downloading at this scale hit NCBI's per-IP rate limit (HTTP 429) under
  8-way parallel downloads; fixed with retry/backoff logic and reduced
  concurrency (4 parallel workers). 179 of 180 genomes downloaded
  successfully — the one failure was a GenBank-only accession without a
  RefSeq pair, which likely wouldn't carry the standardized `gene=` symbols
  the matching approach depends on anyway.
- The core scan logic (`find_divergent_disruptions.py`) was unchanged, just
  refactored so its parsing/matching functions could be imported and reused
  by the Part 4 script below (wrapped the main loop in `run_full_scan()` /
  `if __name__ == "__main__"` instead of bare module-level code).

### Results

| Species | Strain-level hits | Unique gene pairs disrupted | ...with mobile-element keyword signature |
|---|---|---|---|
| *E. coli* | 88 | 18 | 13 |
| *B. subtilis* | 18 | 6 | 2 |
| *B. fragilis* | 1 | 1 | 1 |
| *V. cholerae* | 1 | 1 | 1 |
| *C. jejuni*, *M. tuberculosis*, *P. aeruginosa*, *S. aureus*, *S. pyogenes* | 0 | 0 | 0 |

(counts after the same 800bp–200kb size filter used in Part 1, to exclude
large-scale genome rearrangements rather than discrete insertions)

Going from 4 to 20 genomes per species turned up disruptions in *B. subtilis*
and confirmed *B. fragilis*/*V. cholerae* candidates that the tiny Part 1
sample had missed or only barely caught — exactly the effect a bigger sample
should have. Five species (*C. jejuni*, *M. tuberculosis*, *P. aeruginosa*,
*S. aureus*, *S. pyogenes*) still show zero candidates even at 20 strains;
this is more informative than the Part 1 null result, but still not proof of
absence — see the gene-symbol-matching caveat below.

**Caveat carried over and sharpened:** matching genes across strains by their
RefSeq `gene=` symbol only works for genomes with standardized, RefSeq-style
annotation. Genomes annotated by other pipelines, or with many genes left as
unnamed hypotheticals, are invisible to this method regardless of how many
are sampled. This likely contributes to persistent zero-hit species and is
the main reason a protein-BLAST-based ortholog map (feasible with the
`blastp`/`makeblastdb` already available locally) is the natural next
upgrade, rather than sampling even more genomes with the same method.

## Part 4: Starting from functionally-coupled divergent pairs, then checking for disruption

### Motivation

Part 1/3 start from "any short-spacer divergent pair" and ask whether it's
disrupted; Part 2 showed that almost none of the resulting *E. coli*
candidates are actually functionally coupled the way *torS*/*torT* are. This
part flips the search: start from divergent pairs that are **known or
plausibly functionally coupled**, across bacterial diversity, and ask
whether any of them are ever disrupted by an insertion — directly targeting
a second real example of the *torS*/*torT* phenomenon rather than hoping to
stumble on one.

### Method

Two complementary approaches, both built on `find_functional_pair_disruptions.py`
(imports the Part 1/3 GFF-parsing and circular-span helpers rather than
duplicating them):

1. **Curated list.** A research pass verified (via EcoCyc/RegulonDB/UniProt/
   primary literature) seven specific, precisely-sourced examples of
   divergently-transcribed, functionally-coupled gene pairs — the same
   architecture as *torS*/*torT* (a regulator or accessory factor sharing a
   promoter/operator with the thing it controls): *torS*/*torT*,
   *metE*/*metR* (MetR activates *metE*, methionine biosynthesis),
   *araC*/*araB* (AraC autoregulates and activates *araBAD*),
   *soxR*/*soxS* (SoxR activates *soxS*, oxidative-stress response),
   *acrR*/*acrA* (AcrR represses the *acrAB* efflux operon),
   *mepR*/*mepA* (MepR represses the *mepA* efflux pump, *S. aureus*), and
   *cpxP*/*cpxR* (CpxP divergently transcribed from the *cpxRA* operon it
   helps regulate, conserved into *V. cholerae*). Each pair was checked
   directly (by exact gene symbol) for presence, divergent orientation, and
   disruption across every strain of every species in the dataset — not just
   the species where the relationship was first described, since some of
   these are conserved across many Gammaproteobacteria.
2. **General heuristic.** Rather than hand-curating a pair for every
   species, a broader proxy was applied to the divergent pairs already found
   in Part 1/3: flag any pair where one gene's product annotation names a
   regulator (`TetR|LysR|MarR|AraC|GntR|...family`, "transcriptional
   regulator/repressor/activator", "response regulator", "sigma factor",
   etc.) — the general architecture behind every one of the curated examples
   above. This is a **screen, not proof**: being a regulator doesn't mean it
   regulates *this specific* neighboring gene (Part 2 already demonstrated
   this with *fabR*/*sthA* — FabR is a real repressor, just not of *sthA*).
   So every heuristic hit that showed disruption in ≥1 strain was then
   individually checked against the literature before being reported as a
   real functional link.

### Results

**Curated pairs:**

| Pair | Species checked | Outcome |
|---|---|---|
| *torS*/*torT* | *E. coli* (baseline 82bp) | **Disrupted in 3 of the new random strains** — reconfirms the known HK022-type phenomenon at scale, independent validation that the pipeline recovers the paper's own case |
| *torS*/*torT* | *V. cholerae* (baseline 87bp) | present, not disrupted in this sample |
| *metE*/*metR* | *V. cholerae* (baseline 395bp) | **Disrupted in 1 strain** (GCA_045184175.1): spacer blown up to 14.7kb, intervening genes include an integrase — **a genuine second example of the torS/torT phenomenon**, in a different species, with MetR's real, literature-confirmed regulatory target (*metE*) |
| *metE*/*metR* | *E. coli* (baseline 236bp) | present, not disrupted in this sample |
| *araC*/*araB*, *soxR*/*soxS*, *acrR*/*acrA*, *cpxP*/*cpxR* | *E. coli* | all present with normal short spacers, **not disrupted** in this sample |
| *mepR*/*mepA* | *S. aureus* | not found (likely absent from these particular assemblies or annotated without matching `gene=` symbols — this pair is often plasmid-borne) |

**Heuristic scan + verification:** the heuristic flagged 55–262 regulator-headed
divergent pairs per species (most never disrupted in the sample). Every hit
that *did* show disruption was checked individually:

| Pair | Species | Verdict | Disruption details |
|---|---|---|---|
| *rocD*/*rocR* | *B. subtilis* | **CONFIRMED functional link** — RocR (NtrC/NifA-family, σL-dependent) is the documented activator of the *rocDEF* arginine-degradation operon that *rocD* belongs to (Calogero 1994; Gardan 1995/1997) | Disrupted in 6 of 19 sampled strains: spacer 240bp → ~3.8kb. The inserted genes are *not* a classic phage/IS element — they're the *sdpABC*/*sdpRI* "cannibalism" toxin-immunity-repressor module (Gonzalez-Pastor et al. 2003), present at this exact locus in some *B. subtilis* strains and absent in others. Open question: true horizontal insertion vs. a strain-specific genomic rearrangement placing this native accessory module next to *rocD*/*rocR* — worth follow-up. |
| *yyaT*/*yybA* | *B. subtilis* | **CONFIRMED, pathway-consistent** — comparative-genomics regulon reconstruction (RegPrecise/Rodionov lab) places *yyaT* in YybA's (MarR-family) predicted polyamine-degradation regulon; not yet ChIP/EMSA-validated | Disrupted in 1 of 19 strains: spacer 119bp → 1.5kb, with an **IS4-family transposase** — a clean, compact IS-element insertion |
| *btsS*/*mlrA* | *E. coli* | **NOT CONFIRMED — coincidental neighbor.** MlrA's real, documented target is *csgD* (curli/biofilm regulation, Ogasawara 2010); the primary paper on this exact locus (Behr et al. 2012) explicitly notes mlrA sits upstream of btsS with no functional connection to BtsSR pyruvate sensing | Disrupted in 7 of 19 strains by a recognizable lambdoid prophage (genes named *xisR*, *bet*, *kil*, *gpU/gpG/gpJ*); one strain's inserted phage carries **Shiga toxin genes (*stxA1a*/*stxB1a*)** — a striking, real Stx-converting-phage insertion, just not at a functionally-coupled pair |
| *purT*/*ybfI* | *B. subtilis* | **NOT CONFIRMED** — YbfI's function/regulon is uncharacterized (SubtiWiki/RegPrecise); no evidence it targets *purT* or purine biosynthesis | Disrupted in 5 of 19 strains by a small (~3.6–6.3kb) element with phage holin genes |
| *fabR*/*sthA* | *E. coli* | Already established as a **non-link** in Part 2 (FabR's regulon is *fabA*/*fabB*, not *sthA*) — reappears here only because it also happens to pass the "regulator neighbor" heuristic | (already documented in Part 1/3) |

Also flagged and discarded as circular-genome-coordinate/rearrangement
artifacts (>200kb spans, same issue fixed in Part 1): *pgrR*/*ycjY* and
*nimR*/*nimT* in *E. coli*, *ponA2*/*whiB4* in *M. tuberculosis*.

### Interpretation

Casting a wider net for functionally-coupled divergent pairs paid off:
**two solid new examples of the torS/torT phenomenon** turned up —
*metE*/*metR* in *V. cholerae* (a regulator/target pair, disrupted by a
phage-like integrase-bearing element) and *yyaT*/*yybA* in *B. subtilis* (a
predicted regulon member, disrupted by a clean IS4-family transposase) — plus
a third case, *rocD*/*rocR*, that's a confirmed functional pair but disrupted
by a native accessory toxin-immunity module rather than a classic phage/IS
element, raising an open question about mechanism.

Just as important: the generalized "regulator sitting next to a divergent
partner" heuristic produced two more false leads (*btsS*/*mlrA*,
*purT*/*ybfI*) on top of the one from Part 2 (*fabR*/*sthA*) — a majority
of its "hits" so far have not held up under literature verification. The
heuristic is useful for generating candidates at scale across many species
without hand-curating gene symbols everywhere, but every candidate it
produces needs the same individual verification done here before being
trusted; it cannot substitute for that step.

Taken together with Part 2, the overall picture is: divergent-pair disruption
by mobile/accessory elements is fairly common (108+ candidates across two
scans), but disruption of a pair that is *also* functionally coupled — the
specific phenomenon in the paper — is much rarer, on the order of a handful
of cases even after deliberately searching for it. That rarity is itself an
interesting result: it suggests HK022/*torS*/*torT* isn't a fluke of *E. coli*
biology specifically, but neither is it the default outcome whenever an
element lands between two divergent genes — most such landings appear to be
functionally silent with respect to the two flanking genes' own regulatory
relationship.

## Part 5: How much data would it take to say anything about selection?

### Motivation

Parts 1–4 are descriptive: they establish that insertions into the shared
regulatory regions of divergent pairs happen, and that a handful of those
pairs are genuinely functionally coupled. They do not address the
evolutionary question — is there **selection for or against** elements
landing between divergent, functionally-related genes? This part asks what it
would actually take to answer that, using the base rates we measured rather
than guesses (`analysis/power_analysis.py`, rerunnable; derives its rates
directly from the two result JSONs).

### Empirical base rates (from the 213 usable genomes of Part 3)

| Quantity | Value |
|---|---|
| Divergent pairs monitored (all species) | 4,099 |
| ...regulator-headed (coupling proxy) | 771 |
| ...non-regulator-headed | 3,328 |
| Strain-level disruption hits | 108 |
| **Unique pairs** disrupted in ≥1 strain | **26** (4.2x pseudo-replication collapse) |
| P(disrupted \| any divergent pair) | 0.63% |
| P(disrupted \| regulator-headed) | 0.78% |
| P(disrupted \| non-regulator-headed) | 0.60% |
| P(disrupted \| literature-verified coupled) | 0.31% |

### Headline result: the current data cannot determine even the SIGN of the effect

Running the comparative test two ways, on the *same* genomes, gives opposite
answers depending only on how "functionally related" is defined:

| Definition of "functionally coupled" | Coupled vs non-coupled | Odds ratio (95% CI) | Fisher exact p |
|---|---|---|---|
| Regulator-headed heuristic (unverified) | 1.17% vs 0.51% — **enriched** | 2.30 [1.02, 5.18] | 0.046 |
| Only literature-verified pairs | 0.26% vs 0.51% — **depleted** | 0.51 [0.12, 2.20] | 0.556 |

The unverified heuristic yields a nominally significant *enrichment*; the
verified subset points to *depletion* and is nowhere near significant. This
is the single most important design finding in the whole survey: **the binding
constraint is the functional-coupling annotation, not the number of genomes.**
Scaling to 10,000 genomes while still using the unverified heuristic would buy
a confidently wrong answer with a very small p-value. Part 4 already showed
why: a majority of heuristic hits (*btsS*/*mlrA*, *purT*/*ybfI*,
*fabR*/*sthA*) did not survive literature checking, and those false positives
are exactly what drives the spurious enrichment signal.

### Genomes required (80% power, alpha = 0.05)

| Goal | Genomes needed |
|---|---|
| Tier 1: "does this happen beyond *torS*/*torT*?" | ~180 — **already achieved** (*metE*/*metR*, *yyaT*/*yybA*) |
| Tier 2: rate estimate to ±30% | ~1,000 |
| Tier 2: rate estimate to ±20% | ~2,200 |
| Tier 3: detect 5x depletion (strong purifying selection) | ~700 |
| Tier 3: detect 2x depletion (moderate) | ~2,200 |
| Tier 3: detect 1.5x depletion (weak) | ~5,400 |
| Tier 3: detect 5x enrichment (hotspot/adaptive) | ~130 |

Sampling depth instead of breadth helps: at 100 strains/species (raising the
chance of catching polymorphic, recent insertions) a 2x depletion becomes
detectable with ~1,800 genomes across ~18 species, versus ~2,200 across ~91
species at 20 strains each.

**All of these are floors.** They assume independent insertion events, but our
108 strain-level hits collapse to 26 unique pairs, and strains sharing an
insertion by common descent represent one evolutionary event, not many.
Correcting for phylogeny (counting insertions per independent lineage on a
core-genome tree) plausibly inflates every number above by 2–5x.

**Practical bottom line:** ~500–1,000 genomes for a provisional read if the
effect is strong; ~2,000–5,000 for moderate effects; ~10,000+ after
phylogenetic correction for weak ones. At the download throughput achieved in
Part 3, genome acquisition is 1–2 days of wall-clock — genuinely not the
hard part.

### One provisional conclusion available right now

A 5x effect in either direction would need only ~130–700 genomes to detect,
and we do not see one. So **very strong selection for or against insertion at
functionally-coupled divergent pairs is provisionally excluded**; whatever is
happening is a modest effect, which is precisely why it is expensive to
measure.

### Three design changes that beat brute-force scaling

1. **Replace hand-verification with a scalable coupling annotation** —
   RegPrecise/RegulonDB regulon assignments, STRING interaction scores, or
   KEGG pathway co-membership. Highest-value fix by far; without it the sign
   of the effect stays undetermined at any N.
2. **Flip the denominator.** Instead of "what fraction of coupled pairs get
   disrupted" (a rare-event regime needing huge N), ask "of all insertions
   genome-wide, are they positioned to disrupt shared regulatory elements more
   or less often than expected?" That yields thousands of events from genomes
   already in hand, with spacer length, AT content, and attB-site presence as
   covariates to control for the non-selective drivers of insertion site
   preference.
3. **Use allele frequencies rather than presence/absence.** Insertions at
   coupled pairs that remain polymorphic and never fix within a clade imply
   purifying selection; ones fixed across whole clades imply neutral or
   beneficial. Far more informative per genome, but requires within-species
   population depth (~100 strains each) rather than breadth.

## Part 6: Fixing the annotation layer

### Motivation

Part 5 identified the binding constraint: the sign of the selection effect
flips depending on whether "functionally coupled" is defined by the
regulator-keyword heuristic or by hand literature verification. Hand
verification does not scale past a dozen pairs, and the heuristic is wrong
often enough to invert the answer. So the annotation layer — not genome count
— is what needs fixing first.

### Source choice: STRING, not KEGG

KEGG pathway co-membership was rejected because it captures *metabolic*
co-membership and misses *regulatory* relationships. MetR is a transcriptional
activator, not a member of the methionine-biosynthesis map, so KEGG would
score *metE*/*metR* — our strongest Part 4 finding — as uncoupled. STRING
covers all 9 species, captures regulatory/physical/curated associations, and
yields a **continuous** score, which supports regression rather than a
low-power 2x2.

Downloads are modest: STRING v12 per-species detailed links + aliases for 9
species is **73 MB**, plus 6 MB of proteomes — no bulk database install.

### The circularity trap (most important methodological point)

STRING's headline `combined_score` folds in a **`neighborhood` channel computed
from conserved genomic adjacency**, plus `fusion` and `cooccurence`, which are
likewise genomic-context-derived. Our pairs are adjacent *by construction*, so
scoring them with `combined_score` would reward them for being neighbours and
manufacture the very association we are trying to test.

Channel policy adopted:

| Channels | Use |
|---|---|
| `neighborhood`, `fusion`, `cooccurence` | **excluded always** (genomic-context; circular here) |
| `experimental`, `database`, `textmining` | **primary coupling score** |
| + `coexpression` | sensitivity variant only — not positional, but two genes sharing a divergent promoter are co-regulated *as a consequence of the architecture*, so it is reported separately rather than trusted |

Channels are combined with STRING's own probabilistic formula (remove the
0.041 prior per channel, combine as independent evidence, re-add prior) rather
than a crude sum or max.

A second deliberate choice: pairs that cannot be mapped to STRING are recorded
as **missing (`None`), never as zero**. "No annotation available" and
"annotated as uncoupled" are different states, and conflating them would bias
the selection test toward whichever direction the unmappable pairs happen to lean.

### Problem found: gene-symbol joins fail badly for a third of the species

The first pass joined our genes to STRING by gene symbol. Overall 67% of the
4,099 divergent pairs mapped, but the failures were wildly uneven:

| Species | Pairs mapped by symbol |
|---|---|
| *S. aureus* | 98% |
| *M. tuberculosis*, *P. aeruginosa*, *C. jejuni* | 95% |
| *E. coli* | 84% |
| *B. subtilis* | 77% |
| *S. pyogenes* | **6%** |
| *V. cholerae* | **5%** |
| *B. fragilis* | **1%** |

Diagnosis: STRING's reference strain for those three carries no gene-symbol
aliases — only locus tags, and for *V. cholerae* 2000-era GenBank protein
accessions (`AAF93179.1`) that share no namespace with our RefSeq `WP_`
accessions. STRING's `preferred_name` field was checked as a cheap fallback and
supplies real symbols for only 10–30% of proteins in those species. This is
precisely the standardized-symbol fragility already flagged as a limitation in
Part 3 — and it matters, because *V. cholerae* holds our best finding
(*metE*/*metR*).

### Fix: sequence-based ortholog mapping

Replaced the name join with homology: for each species, the proteins of
genes appearing in a divergent pair are searched with `blastp` against that
species' STRING proteome, taking the best hit by bitscore under strict
thresholds (≥70% identity, ≥70% query coverage, E ≤ 1e-20). One-directional
best-hit is defensible here because query and subject are the *same species*,
so the true ortholog is a near-identical, unambiguous top hit; reciprocal best
hits would be the stricter choice for cross-species work. The annotation
script prefers the ortholog map and falls back to symbol, then locus tag.

Beyond rescuing the three failing species, this **removes the gene-symbol
dependency** that has limited every part of the survey since Part 1.

Because this is ~CPU-hours of BLAST, it is submitted to an SGE compute node
(`feder-short.q`, 8 slots on one host via the `serial` PE, threads bound to
`$NSLOTS`) rather than run on the shared login node —
see `analysis/submit_ortholog_map.sh`.

### Ortholog mapping result

The BLAST job (30 s wall-clock on 8 slots — my "CPU-hours" estimate was far too
pessimistic, since only ~1,000 query proteins per species are needed) rescued
all three failing species and improved every other one:

| Species | Symbol join | BLAST orthologs |
|---|---|---|
| *B. fragilis* | 1% | **88%** |
| *V. cholerae* | 5% | **88%** |
| *S. pyogenes* | 6% | **87%** |
| *E. coli* | 84% | 90% |
| *B. subtilis* | 77% | 95% |
| *M. tuberculosis*, *P. aeruginosa*, *C. jejuni* | 95% | 97% |
| *S. aureus* | 98% | 94%* |

Overall divergent-pair coverage rose from **67% → 89%** (3,666 of 4,099), with
no species below 78%. (*For *S. aureus*, a few genes resolve via locus-tag
fallback rather than BLAST, so the pair-level rate shifts slightly.)

### Validation against the hand-verified labels — PASSED

Parts 2 and 4 produced **19 usable labels** (9 verified coupled, 10 verified
*not* coupled — Part 2 checked nine pairs and found no functional link for
any). Scoring them with the new annotation:

| | n | median coupling score |
|---|---|---|
| verified **coupled** | 9 | **897** |
| verified **not** coupled | 10 | **45** |

**AUC = 0.872** — i.e. a truly-coupled pair outranks a truly-uncoupled one ~87%
of the time (Mann-Whitney over all 9x10 comparisons; threshold-free, so it
measures ranking quality rather than accuracy at some cutoff). With n=19 the
interval around 0.872 is wide, so read it as "clearly informative," not
"precisely 0.87."

**Two informative failures**, both false negatives, and both in the same
direction:

- *V. cholerae* **metE/metR scored only 69** despite being a textbook
  regulator/target pair — STRING's *V. cholerae* annotation is sparse.
- *B. subtilis* **yyaT/yybA scored 41** (the prior-only floor) — its regulon
  assignment is a comparative-genomics *prediction* (RegPrecise) that STRING
  does not ingest.

Both are consequences of a real weakness: **the score is carried almost
entirely by the `textmining` channel** (`experimental` and `database` are 0 for
nearly every pair here), so coupling detection tracks *literature depth*. That
biases toward under-calling coupling in non-model organisms — which is exactly
where our two novel findings live. Any future use of this layer should treat
its negatives in understudied species as unreliable.

On the excluded channels: `neighborhood` alone gives AUC 0.656 on the same
labels. That is well below the coupling score's 0.872 but not the ~0.5 of pure
noise — conserved adjacency genuinely correlates with functional coupling
(that is why STRING ships the channel). It remains the wrong channel *here*,
because every pair in this study is adjacent in the reference genome by
construction, so including it would partly score pairs on the property that
defines the study population.

### Payoff: the selection test no longer flips

Re-running the Part 5 comparison on the 3,666 mapped pairs, calling a pair
coupled at STRING's conventional medium-confidence cutoff (score ≥ 400):

| Coupling definition | Disruption rate | Odds ratio (95% CI) | Fisher p |
|---|---|---|---|
| Regulator-keyword heuristic | 1.17% vs 0.51% | 2.30 [1.02, 5.18] | **0.046** |
| Hand-verified pairs only | 0.26% vs 0.51% | 0.51 [0.12, 2.20] | 0.556 |
| **STRING score, validated** | **0.80% vs 0.70%** | **1.15 [0.39, 3.35]** | **0.773** |

The nominally significant enrichment was an **annotation artifact**. With a
validated coupling score there is no detectable effect in either direction, and
the confidence interval [0.39, 3.35] is exactly what Part 5's power analysis
predicted for this sample size: only effects stronger than roughly 3x are
excluded.

Two useful side results:

- **Spacer length is not a confound** — median 174 bp (coupled) vs 180 bp
  (uncoupled), so the obvious non-selective driver of insertion is balanced
  across the comparison groups.
- Disruptions remain overwhelmingly concentrated in *E. coli* (18 of 26) and
  *B. subtilis* (6), so the test is effectively powered by two species.

### What this part did and did not settle

**Settled:** the annotation layer is fixed and, importantly, *validated* rather
than merely automated. The Part 5 sign ambiguity is resolved — the honest
estimate is "no effect detectable," not "enriched."

**Not settled:** the sample is still far too small (26 independent events). The
annotation fix makes the *estimate* trustworthy; it does not make the *sample*
adequate. Part 5's genome targets still stand.

## Files

- `analysis/fetch_genomes.py` — pulls genome assemblies from NCBI Datasets API
  (as of Part 3: random sample of 20 strains/species from a pool of up to
  1000 complete genomes, fixed random seed)
- `analysis/find_divergent_disruptions.py` — core scan (divergent-pair
  discovery, circular-aware cross-strain comparison, mobile-element keyword
  check); parsing/matching functions importable via `run_full_scan()` guard
- `analysis/find_functional_pair_disruptions.py` — Part 4 script: curated
  functionally-coupled pairs + general regulator-heuristic scan, reusing the
  above module
- `analysis/power_analysis.py` — Part 5 design/power analysis; derives base
  rates from the result JSONs and reports genomes needed per goal
- `analysis/annotate_coupling_string.py` — Part 6 annotation layer: attaches
  continuous STRING coupling scores to every divergent pair, excluding the
  circular genomic-context channels
- `analysis/build_ortholog_map.py` — Part 6 BLAST ortholog mapper (replaces the
  fragile gene-symbol join); run via the SGE submit script below
- `analysis/submit_ortholog_map.sh` — SGE submit script (`feder-short.q`, 8
  slots); keeps the BLAST work off the login node
- `data/string/` — STRING v12 per-species links/aliases/proteomes (73 MB + 6 MB)
- `data/refprot/` — reference-genome protein FASTAs (BLAST queries)
- `data/coupling_scores.json` / `.tsv` — per-pair coupling scores (symbol-join
  pass; to be regenerated with the ortholog map)
- `analysis/validate_coupling.py` — scores the 19 hand-verified labels, reports
  AUC, and quantifies the `neighborhood` circularity
- `analysis/test_selection.py` — re-runs the selection test on validated scores
- `data/ortholog_map.json` — BLAST ortholog map (89% pair coverage)
- `DECK.html` — 2-3 min lightning talk (self-contained; arrow keys, T for theme)
- `logs/` — SGE job logs
- `data/species_list.txt`, `data/genome_manifest.tsv` — exact species/accessions used
- `data/genomes/` — downloaded GFF3 + FASTA per genome (180 requested / 179
  retrieved / 213 present including Part 1 holdovers, as of Part 3).
  **Gitignored** — large and fully re-fetchable via `fetch_genomes.py`
  plus `genome_manifest.tsv`
- `data/divergent_disruption_results.json` — full raw output of the Part 1/3 scan
- `data/functional_pair_results.json` — full raw output of the Part 4 scan
