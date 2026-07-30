# Homology search for phage-modulated divergent-promoter architectures

## Reference

Carey et al. 2019, "Phage integration alters the respiratory strategy of its host"
(`Carey et al. - 2019 - Phage integration alters the respiratory strategy .pdf`).

- Supp 1 (`Carey et al. - 2019 Supplement 1.docx`): table of 46 *E. coli*
  strains carrying a prophage in the *torT*–*torS* intergenic region, with
  phage identity, phage accession, and expanded intergenic distance.
- Supp 2 (`Carey et al. - 2019 Supplement 2.pdf`): 220-bp nucleotide MSA of the
  intergenic region across all 46 strains + MG1655::HK022, annotated with
  `attL_HK022` (B, O, P′), the phage-encoded TSS, and the `int` gene position.

The system under study: two genes (*torT*, *torS*) coregulated from an
overlapping/back-to-back intergenic promoter region; phage integration
between them expands that region and rewires expression.

## Data on hand

### MMseqs2 clustered protein DBs
Two versions exist under `/fh/working/srivatsan_s/databases/mmseqs2/`,
both amino-acid (`.dbtype = 0`) with the same header scheme
`omg_row<N>_cds<M>` but indexing different Arrow splits:

**Full — `mmseqs2/omg_cds/`** (~3.9 TB on disk)
- Source `omg_cds.fasta`: **996 GB** (3.28 billion residues); built from the
  full `train/` split (3833 shards, ~all OMG contigs). Verified:
  `omg_row0_cds0` = `MALTKVEKRNRIKRR...` = `train[0]['CDS_seqs'][0]`.
- `omg_cds_clustered90_seqs*`: **171 GB** representative seqs at 90% ID
  (~10× the 10M subset).
- `omg_cds_db_gpu*`: GPU-padded copy of the **un-clustered** full DB
  (~924 GB). **No GPU-padded copy of the cluster representatives exists.**
  So on GPU you'd search the un-clustered set (losing dedup); on CPU you
  can search the cluster reps.

**10M subset — `mmseqs2/omg_cds_10M/`** (~186 GB on disk;
symlinked from `~/working/databases/mmseqs2/omg_cds_10M/`)
- Source `omg_cds.fasta`: 34 GB, built from the `train_10M_random/` split
  (200 shards, ~10M contigs, random sample). Verified: `omg_row0_cds0`
  = `MFYKIIA...` = `train_10M_random[0]['CDS_seqs'][0]`.
- `omg_cds_clustered90_seqs*`: 18 GB cluster representatives.
- `omg_cds_clustered90_seqs_pad*`: 18 GB GPU-padded copy **of the cluster
  representatives** — ideal for a fast pilot with `mmseqs search --gpu 1`.

**Recommended order:** pilot on the 10M subset (GPU search on padded cluster
reps is fast, taxonomic breadth is already large), then rerun on the full
DB after the pipeline is validated. Note: the full DB uses different row
indexing (into `train/`), so IGS lookups (step 5) must switch splits.

- MMseqs2 available on the cluster via `module load MMseqs2/13-45111-gompi-2021b`.

### Underlying OMG Arrow dataset (unlocks promoter DNA)
Path: `/fh/working/srivatsan_s/databases/data/OMG/` — HuggingFace-style Arrow
dataset with two splits:

- `train_10M_random/` — 200 shards, ~10M contigs (random sample). Source
  of the `omg_cds_10M/` mmseqs DB. Header `omg_row<N>_cds<M>` in that DB
  maps to `train_10M_random[N]['CDS_seqs'][M]`.
- `train/` — 3833 shards, full (~29× larger) unsubsampled OMG. Source of
  the full `omg_cds/` mmseqs DB. Header `omg_row<N>_cds<M>` in that DB
  maps to `train[N]['CDS_seqs'][M]`.

Schema per row (contig):

| field | type | contents |
|---|---|---|
| `CDS_position_ids` | list<int32> | position index of each CDS (odd numbers) |
| `IGS_position_ids` | list<int32> | position index of each intergenic segment (even numbers) |
| `CDS_ids` | list<string> | `{taxon}\|{contig}\|CDS\|{gene}\|{strand}\|{start}:{end}` |
| `IGS_ids` | list<string> | `{taxon}\|{contig}\|IG\|{ig_id}\|{strand}\|{start}:{end}` |
| `CDS_seqs` | list<large_string> | **protein** sequences |
| `IGS_seqs` | list<large_string> | **nucleotide** sequences (intergenic DNA) |
| `CDS_orientations` | list<bool> | per-CDS strand |

Positions interleave as `IGS0, CDS0, IGS1, CDS1, IGS2, CDS2, ...`, so the
IGS *between* CDS_i and CDS_{i+1} is at IGS index `i+1`. That means once
we identify a candidate torT/torS homolog pair by mmseqs hit → contig row +
CDS indices, we can pull the promoter DNA (with original genomic
coordinates and taxon/contig IDs) with a direct Arrow lookup — no external
NCBI fetch needed. **This unblocks promoter-level analysis.**

### Local python env
`micromamba activate omg_search` (created 2026-07-30) — includes
`pyarrow`, `datasets`, `biopython`, `pandas`, `duckdb`, `tqdm`.
Note: a stale `~/.local/lib/python3.11/site-packages/pandas` shadows the
conda pandas; either delete that dir or run with `PYTHONNOUSERSITE=1`.

## What is still missing

1. **Query proteins.** Neither supplement contains protein FASTA. Pull from
   NCBI:
   - *torT* (periplasmic TMAO-binding, ~338 aa) and *torS* (hybrid sensor
     kinase, ~914 aa) from *E. coli* MG1655 (U00096 / P36664 / P39453).
   - Integrase (`int`) protein for each phage in Supp 1: Sf101
     (NC_027398), BP-4795 (NC_004813), CDT-1Φ (NC_009514), DE3 (NC_042057),
     GF-2 (NC_026611), HK620 (NC_002730), SEN34 (NC_028699), Sf6
     (NC_005344), YYZ-2008 (NC_011356), λ (NC_001416), HK022 (NC_002166).
2. **Phage-gene profile set (optional, later).** Beyond just integrase, a
   small HMM/profile library of core phage genes (terminase, capsid, tail,
   portal) would tighten the "prophage inserted between the pair" call.
   Recommend PHROGs; not currently local.

Note: the earlier "we don't have promoter DNA" gap is *closed* — see IGS
fields above.

## Plan

### 0. Prerequisites

- `mkdir -p queries/ results/ tmp/ logs/ scripts/`
- On rhino: `ml MMseqs2/13-45111-gompi-2021b && mmseqs version`
- On rhino: `micromamba activate omg_search` (run with `PYTHONNOUSERSITE=1`
  until the stale user-site pandas is cleaned).
- GPU node for search: `sbatch --gres=gpu:1` (gizmo k-gen RTX 2080ti 11 GB
  is enough for the 18 GB padded 10M-cluster DB via disk-backed access;
  gizmo-harmony L40S 44 GB is the target for the full DB later).
  Per-account cap: 20 concurrent GPUs.

### 1. Assemble query FASTAs

Pilot design (agreed 2026-07-30):

- `queries/torT_torS.faa` — MG1655 TorT (P36664) and TorS (P39453) only.
  Ortholog broadening deferred.
- `queries/phage_integrases.faa` — one Int protein per phage in Supp 1
  (11 sequences: HK022, λ, Sf6, Sf101, BP-4795, CDT-1Φ, DE3, GF-2,
  HK620, SEN34, YYZ-2008). No PHROG HMM expansion for pilot.

### 2. Build query DBs and run search on the clustered OMG DB

```
mmseqs createdb queries/torT_torS.faa queries/torT_torS_db
mmseqs createdb queries/phage_integrases.faa queries/phage_int_db

mmseqs search queries/torT_torS_db \
              ~/working/databases/mmseqs2/omg_cds_10M/omg_cds_clustered90_seqs_pad \
              results/torT_torS_vs_omg  tmp/ \
              -s 7.5 -e 1e-10 --max-seqs 5000 --gpu 1

mmseqs search queries/phage_int_db \
              ~/working/databases/mmseqs2/omg_cds_10M/omg_cds_clustered90_seqs_pad \
              results/phage_int_vs_omg  tmp/ \
              -s 7.5 -e 1e-10 --max-seqs 20000 --gpu 1
```

- Convert to `.m8` TSV with `mmseqs convertalis`, keeping the target header
  field (`omg_row<N>_cds<M>`) — that's the key for step 3.
- Sensitivity `-s 7.5` because query/target divergence across metagenomes is
  high; `--gpu 1` uses the `_pad` DB. Drop `--gpu 1` if no GPU node available.

### 3. Reconstruct synteny from cluster hits (Arrow-backed)

For each contig row `N`:

- Parse target headers → `(N, M)` pairs for each query class (torT, torS,
  integrase).
- Test the paired-and-adjacent condition:
  - both torT-like and torS-like on the same row `N`
  - `|M_torT - M_torS| <= K` (start with `K = 3`, also record strictly
    adjacent `K = 1`)
- For each candidate pair, list all CDS between them and check whether any
  are in the integrase (or later, phage-gene) hit set.

Four architecture classes per candidate locus:

1. **A / baseline** — adjacent torT/torS pair, no intervening phage hit.
2. **B / interrupted** — torT/torS pair separated by 1–N CDS with **≥1**
   integrase (or phage-gene) hit between them. The Carey analog.
3. **C / diverged** — torT/torS pair separated by many CDS with no phage
   signature (either genuine rearrangement or assembly break).
4. **B′ / half-pair + phage** — contig has torT (or torS) adjacent to a
   phage-integrase hit but the second anchor is missing from the contig
   (likely truncated at the contig edge). Counted as a candidate per
   decision 2026-07-30, but flagged separately for review.

### 4. Rank, deduplicate, and store

- Load `train_10M_random` shards with pyarrow (map row `N` → the target
  shard by `N // rows_per_shard`; confirm rows/shard on the first pass).
- Expand mmseqs cluster representatives to members with
  `mmseqs createtsv omg_cds_clustered90 omg_cds_clustered90 members.tsv` so
  we don't undercount source contigs.
- Score each locus on: bidirectional best-hit for torT and torS,
  integrase bit-score, orientation (Carey pair is divergent — CDS strands
  should be `-, +` for `torT, torS`), and phylogenetic breadth (count of
  distinct taxa via the leading `{taxon}` field of `CDS_ids`).
- Persist to DuckDB (`results/candidates.duckdb`) with tables
  `mmseqs_hits`, `pair_candidates`, `pair_igs` (see step 5). DuckDB gives
  cheap joins between the two searches and the Arrow rows.

### 5. Promoter-level follow-up (Arrow lookup — now feasible)

For each candidate locus from step 4:

- Look up `train_10M_random[N]` → grab `CDS_ids`, `CDS_orientations`,
  `IGS_ids`, `IGS_seqs`, and position IDs.
- Recover the promoter DNA: for a candidate pair at CDS indices
  `i, j` (with `i < j`), the intergenic sequences at IGS indices
  `i+1 .. j` are the phage-modulated region. For the un-integrated class
  (`j = i+1`), a single IGS is the shared divergent promoter.
- Compare to the Supp 2 alignment: is there a conserved overlapping
  divergent-promoter block? Does an `attB`-like core sit in it? Does the
  interrupted class show an inserted `attL/attR` signature flanking the
  phage-gene stretch?
- MSA the un-integrated-class IGS with MAFFT / MUSCLE; motif-scan with
  MEME/FIMO for divergent −10/−35 pairs on opposite strands.

### 6. Scale up and generalize (stretch)

- **Full DB rerun.** After the pipeline is validated on the 10M subset,
  rerun steps 2–5 against `mmseqs2/omg_cds/` (~29× more contigs, ~10× more
  cluster representatives). Two operational choices:
  - CPU on `omg_cds_clustered90_seqs` (dedup preserved, slower).
  - GPU on `omg_cds_db_gpu` (fast but searches the un-clustered 924 GB
    set — expect much larger hit lists and heavier post-processing).
  Step 5 IGS lookup must switch to the `train/` split (row indexing
  differs from the 10M subset).
- **Other divergent-promoter pairs.** The Carey mechanism is not specific
  to torT/torS. Swap in other overlapping-promoter pairs (e.g. *lexA/dinB*,
  *ompR/envZ*, two-component sensor/regulator pairs) and rerun. That turns
  the pipeline into a general screen for "phage-interruptible coregulated
  pairs" across metagenomes.

## Decisions log (2026-07-30)

- Query set: MG1655 torT + torS only for the pilot; ortholog broadening
  and PHROG integrase HMM deferred.
- Compute: GPU search on the 10M padded cluster DB via Slurm.
- Scope: pilot on `omg_cds_10M/` first, then rerun on full `omg_cds/`.
- Class B′ (half-pair + adjacent phage) counted as candidate but flagged.
- Promoter definition: use IGS boundaries as OMG defines them (from CDS
  calls); acknowledge that this is CDS-call-dependent.
