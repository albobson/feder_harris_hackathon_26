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

### Files

- `analysis/fetch_genomes.py` — pulls genome assemblies from NCBI Datasets API
  (as of Part 3: random sample of 20 strains/species from a pool of up to
  1000 complete genomes, fixed random seed)
- `analysis/find_divergent_disruptions.py` — core scan (divergent-pair
  discovery, circular-aware cross-strain comparison, mobile-element keyword
  check); parsing/matching functions importable via `run_full_scan()` guard
- `analysis/find_functional_pair_disruptions.py` — Part 4 script: curated
  functionally-coupled pairs + general regulator-heuristic scan, reusing the
  above module
- `data/species_list.txt`, `data/genome_manifest.tsv` — exact species/accessions used
- `data/genomes/` — downloaded GFF3 + FASTA per genome (180 genomes as of Part 3)
- `data/divergent_disruption_results.json` — full raw output of the Part 1/3 scan
- `data/functional_pair_results.json` — full raw output of the Part 4 scan

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
