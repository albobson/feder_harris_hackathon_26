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
- `analysis/find_divergent_disruptions.py` — core scan (divergent-pair
  discovery, circular-aware cross-strain comparison, mobile-element keyword
  check)
- `data/species_list.txt`, `data/genome_manifest.tsv` — exact species/accessions used
- `data/genomes/` — downloaded GFF3 + FASTA per genome
- `data/divergent_disruption_results.json` — full raw output of the scan
