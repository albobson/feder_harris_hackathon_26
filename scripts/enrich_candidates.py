"""Step 4-5 of the plan: enrich candidate loci with taxonomy/orientation
and pull the intergenic (promoter) DNA out of the OMG Arrow shards.

Input:  results/candidate_loci.tsv     (from classify_candidates.py)
        results/half_pair_loci.tsv     (from classify_candidates.py)
        /fh/working/srivatsan_s/databases/data/OMG/train_10M_random/*.arrow
Output: results/candidate_loci_enriched.tsv
        results/half_pair_loci_enriched.tsv
        results/class_A_igs.fasta          (divergent-promoter DNA per class-A pair)
        results/class_B_igs.fasta          (phage-spanning IGS per class-B pair, if any)
        cache/omg_shard_index.json         (row -> shard cache, ~5 s to rebuild)

Contig row N in the mmseqs headers `omg_row<N>_cds<M>` maps to
`train_10M_random[N]`. Each row has CDS_position_ids and IGS_position_ids
which together tile 0..N-1 across both feature types. IGS is stored ONLY
when there's an actual intergenic gap between adjacent CDS calls — abutting
or overlapping ORFs get no IGS entry, so `n_IGS != n_CDS + 1` in general.

To find IGS between CDS[i] and CDS[j] (i < j): look up
CDS_position_ids[i]=P_i and CDS_position_ids[j]=P_j, then collect IGS
indices k where P_i < IGS_position_ids[k] < P_j.

For class A (adjacent torT/torS, gap=1), there is at most one IGS between —
the shared divergent-promoter DNA. For class B (phage-interrupted), multiple
IGS segments span the insertion.

Run:
    PYTHONNOUSERSITE=1 ~/micromamba/envs/omg_search/bin/python \\
        scripts/enrich_candidates.py
"""

from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import pyarrow as pa
import pyarrow.ipc as ipc

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
CACHE = REPO / "cache"
SHARD_DIR = Path("/fh/working/srivatsan_s/databases/data/OMG/train_10M_random")


def build_shard_index() -> list[tuple[int, int, Path]]:
    """Return [(cum_start, cum_end_exclusive, path), ...] sorted by cum_start.
    Cached to cache/omg_shard_index.json so subsequent runs are cheap."""
    CACHE.mkdir(exist_ok=True)
    cache = CACHE / "omg_shard_index.json"
    if cache.exists():
        entries = json.loads(cache.read_text())
        return [(s, e, Path(p)) for s, e, p in entries]

    shards = sorted(p for p in SHARD_DIR.iterdir()
                    if p.name.startswith("data-") and p.suffix == ".arrow")
    idx: list[tuple[int, int, Path]] = []
    cum = 0
    for shard in shards:
        with pa.memory_map(str(shard), "r") as src:
            reader = ipc.open_stream(src)
            nrows = sum(b.num_rows for b in reader)
        idx.append((cum, cum + nrows, shard))
        cum += nrows
    print(f"built shard index: {len(idx)} shards, {cum} rows total",
          file=sys.stderr)
    cache.write_text(json.dumps([(s, e, str(p)) for s, e, p in idx]))
    return idx


def shard_for_row(idx: list[tuple[int, int, Path]], row: int) -> tuple[Path, int]:
    """Binary search: return (shard_path, local_offset)."""
    lo, hi = 0, len(idx)
    while lo < hi:
        mid = (lo + hi) // 2
        s, e, _ = idx[mid]
        if row < s:
            hi = mid
        elif row >= e:
            lo = mid + 1
        else:
            return idx[mid][2], row - s
    raise IndexError(f"row {row} out of range")


def load_row(shard_path: Path, offset: int) -> dict:
    """Return the row as a python dict (single-row slice)."""
    with pa.memory_map(str(shard_path), "r") as src:
        reader = ipc.open_stream(src)
        seen = 0
        for batch in reader:
            n = batch.num_rows
            if offset < seen + n:
                local = offset - seen
                tbl = pa.Table.from_batches([batch]).slice(local, 1)
                row = {k: v[0] for k, v in tbl.to_pydict().items()}
                return row
            seen += n
    raise IndexError(f"offset {offset} not found in {shard_path}")


def taxon_from_cds_id(cds_id: str) -> str:
    """`7000000126|C1821366|CDS|...` -> `7000000126`."""
    return cds_id.split("|", 1)[0]


def contig_from_cds_id(cds_id: str) -> str:
    """`7000000126|C1821366|CDS|...` -> `C1821366`."""
    parts = cds_id.split("|")
    return parts[1] if len(parts) > 1 else ""


def genomic_coords(cds_id: str) -> tuple[int, int] | None:
    """`...|+|84:437` -> (84, 437). Returns None if unparseable."""
    tail = cds_id.rsplit("|", 1)[-1]
    if ":" not in tail:
        return None
    try:
        a, b = tail.split(":")
        return int(a), int(b)
    except ValueError:
        return None


def igs_between(row: dict, cds_i: int, cds_j: int) -> list[int]:
    """Return IGS indices whose position_id sits strictly between
    CDS_position_ids[cds_i] and CDS_position_ids[cds_j]."""
    cds_pos = row["CDS_position_ids"]
    igs_pos = row["IGS_position_ids"]
    lo, hi = sorted((cds_pos[cds_i], cds_pos[cds_j]))
    return [k for k, p in enumerate(igs_pos) if lo < p < hi]


def enrich_paired(idx, paired_rows: list[dict]) -> tuple[list[dict], list[dict], list[dict]]:
    """Load Arrow rows for each candidate and return
    (enriched_rows, class_A_fasta_records, class_B_fasta_records)."""
    enriched: list[dict] = []
    fasta_A: list[dict] = []
    fasta_B: list[dict] = []

    # Group candidates by shard so we touch each shard at most once.
    by_shard: dict[Path, list[dict]] = defaultdict(list)
    for cand in paired_rows:
        row = int(cand["omg_row"])
        try:
            shard, local = shard_for_row(idx, row)
        except IndexError:
            cand["_error"] = f"row {row} out of shard range"
            enriched.append(cand)
            continue
        cand["_shard"] = shard
        cand["_local"] = local
        by_shard[shard].append(cand)

    for shard, cands in by_shard.items():
        cands_sorted = sorted(cands, key=lambda c: c["_local"])
        with pa.memory_map(str(shard), "r") as src:
            reader = ipc.open_stream(src)
            batch_start = 0
            pending = cands_sorted[:]
            for batch in reader:
                batch_end = batch_start + batch.num_rows
                take = [c for c in pending if batch_start <= c["_local"] < batch_end]
                if take:
                    tbl = pa.Table.from_batches([batch])
                    for cand in take:
                        local_in_batch = cand["_local"] - batch_start
                        row = {k: v[0] for k, v in tbl.slice(local_in_batch, 1).to_pydict().items()}
                        _enrich_one(cand, row, enriched, fasta_A, fasta_B)
                pending = [c for c in pending if c["_local"] >= batch_end]
                if not pending:
                    break
                batch_start = batch_end

    return enriched, fasta_A, fasta_B


def _enrich_one(cand: dict, row: dict,
                enriched: list, fasta_A: list, fasta_B: list) -> None:
    tort_i = int(cand["torT_cds"])
    tors_i = int(cand["torS_cds"])
    lo, hi = min(tort_i, tors_i), max(tort_i, tors_i)
    orientations = row["CDS_orientations"]
    cds_ids = row["CDS_ids"]

    if tort_i >= len(cds_ids) or tors_i >= len(cds_ids):
        cand["_error"] = "CDS index out of range for this row"
        enriched.append(cand)
        return

    tort_strand = "+" if orientations[tort_i] else "-"
    tors_strand = "+" if orientations[tors_i] else "-"
    divergent = orientations[tort_i] != orientations[tors_i]

    taxon = taxon_from_cds_id(cds_ids[tort_i])
    contig = contig_from_cds_id(cds_ids[tort_i])

    # Genomic gap between the two CDS calls (independent sanity check)
    tort_coord = genomic_coords(cds_ids[tort_i])
    tors_coord = genomic_coords(cds_ids[tors_i])
    genomic_gap = ""
    if tort_coord and tors_coord:
        # gap between end of upstream and start of downstream in genome
        tort_a, tort_b = tort_coord
        tors_a, tors_b = tors_coord
        if tort_a < tors_a:
            genomic_gap = tors_a - tort_b
        else:
            genomic_gap = tort_a - tors_b

    cand.update({
        "taxon": taxon,
        "contig": contig,
        "n_cds_in_contig": len(cds_ids),
        "torT_strand": tort_strand,
        "torS_strand": tors_strand,
        "divergent": divergent,
        "torT_cds_id": cds_ids[tort_i],
        "torS_cds_id": cds_ids[tors_i],
        "genomic_gap_bp": genomic_gap,
    })

    # IGS via position-ID lookup (positions interleave CDS and IGS on 0..N-1)
    igs_seqs = row["IGS_seqs"]
    ig_indices = igs_between(row, tort_i, tors_i)
    total_igs_bp = sum(len(igs_seqs[k]) for k in ig_indices)
    cand["n_igs_between"] = len(ig_indices)
    cand["total_igs_bp"] = total_igs_bp

    if cand["class"] == "A" and ig_indices:
        # single divergent-promoter IGS
        k = ig_indices[0]
        seq = igs_seqs[k]
        header = (f"omg_row{cand['omg_row']}_igs{k}"
                  f" taxon={taxon} contig={contig}"
                  f" torT_cds={tort_i}({tort_strand})"
                  f" torS_cds={tors_i}({tors_strand})"
                  f" len={len(seq)}"
                  f" divergent={divergent}")
        fasta_A.append({"header": header, "seq": seq})
    elif cand["class"] == "B" and ig_indices:
        for k in ig_indices:
            seq = igs_seqs[k]
            header = (f"omg_row{cand['omg_row']}_igs{k}"
                      f" taxon={taxon} contig={contig}"
                      f" between torT_cds={tort_i} and torS_cds={tors_i}"
                      f" len={len(seq)}")
            fasta_B.append({"header": header, "seq": seq})

    enriched.append(cand)


def enrich_half_pair(idx, hp_rows: list[dict]) -> list[dict]:
    """Enrich half-pair (B') records with taxon and orientation of the anchor."""
    by_shard: dict[Path, list[dict]] = defaultdict(list)
    enriched: list[dict] = []
    for cand in hp_rows:
        row = int(cand["omg_row"])
        try:
            shard, local = shard_for_row(idx, row)
        except IndexError:
            cand["_error"] = f"row {row} out of shard range"
            enriched.append(cand)
            continue
        cand["_shard"] = shard
        cand["_local"] = local
        by_shard[shard].append(cand)

    for shard, cands in by_shard.items():
        cands_sorted = sorted(cands, key=lambda c: c["_local"])
        with pa.memory_map(str(shard), "r") as src:
            reader = ipc.open_stream(src)
            batch_start = 0
            pending = cands_sorted[:]
            for batch in reader:
                batch_end = batch_start + batch.num_rows
                take = [c for c in pending if batch_start <= c["_local"] < batch_end]
                if take:
                    tbl = pa.Table.from_batches([batch])
                    for cand in take:
                        local_in_batch = cand["_local"] - batch_start
                        row = {k: v[0] for k, v in tbl.slice(local_in_batch, 1).to_pydict().items()}
                        anchor_i = int(cand["anchor_cds"])
                        int_i = int(cand["integrase_cds"])
                        cds_ids = row["CDS_ids"]
                        orientations = row["CDS_orientations"]
                        if anchor_i >= len(cds_ids) or int_i >= len(cds_ids):
                            cand["_error"] = "CDS index out of range"
                            enriched.append(cand)
                            continue
                        cand.update({
                            "taxon": taxon_from_cds_id(cds_ids[anchor_i]),
                            "contig": contig_from_cds_id(cds_ids[anchor_i]),
                            "n_cds_in_contig": len(cds_ids),
                            "anchor_strand": "+" if orientations[anchor_i] else "-",
                            "integrase_strand": "+" if orientations[int_i] else "-",
                            "anchor_cds_id": cds_ids[anchor_i],
                            "integrase_cds_id": cds_ids[int_i],
                        })
                        enriched.append(cand)
                pending = [c for c in pending if c["_local"] >= batch_end]
                if not pending:
                    break
                batch_start = batch_end

    return enriched


def write_tsv(rows: list[dict], path: Path, base_fields: list[str],
              add_fields: list[str]) -> None:
    fields = base_fields + [f for f in add_fields if f not in base_fields]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_fasta(records: list[dict], path: Path) -> None:
    with path.open("w") as fh:
        for r in records:
            fh.write(f">{r['header']}\n")
            seq = r["seq"]
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def main() -> int:
    idx = build_shard_index()

    with (RESULTS / "candidate_loci.tsv").open() as fh:
        paired = list(csv.DictReader(fh, delimiter="\t"))
    with (RESULTS / "half_pair_loci.tsv").open() as fh:
        hp = list(csv.DictReader(fh, delimiter="\t"))
    print(f"loaded {len(paired)} paired candidates and {len(hp)} half-pair records",
          file=sys.stderr)

    enriched_paired, fasta_A, fasta_B = enrich_paired(idx, paired)
    enriched_hp = enrich_half_pair(idx, hp)

    base_paired = ["omg_row", "class", "taxon", "contig", "n_cds_in_contig",
                   "torT_cds", "torS_cds", "gap", "torT_strand", "torS_strand",
                   "divergent", "torT_bits", "torT_pident", "torS_bits",
                   "torS_pident", "n_integrase_between",
                   "integrase_families_between", "integrase_cds_between",
                   "n_igs_between", "total_igs_bp", "genomic_gap_bp",
                   "torT_cds_id", "torS_cds_id", "notes"]
    base_hp = ["omg_row", "anchor", "taxon", "contig", "n_cds_in_contig",
               "anchor_cds", "integrase_cds", "gap_to_integrase",
               "anchor_strand", "integrase_strand", "integrase_family",
               "anchor_bits", "anchor_pident", "integrase_bits",
               "anchor_cds_id", "integrase_cds_id"]

    write_tsv(enriched_paired, RESULTS / "candidate_loci_enriched.tsv",
              base_paired, [])
    write_tsv(enriched_hp, RESULTS / "half_pair_loci_enriched.tsv", base_hp, [])
    write_fasta(fasta_A, RESULTS / "class_A_igs.fasta")
    write_fasta(fasta_B, RESULTS / "class_B_igs.fasta")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print("\n=== paired candidates by class ===", file=sys.stderr)
    class_ct: dict[str, int] = defaultdict(int)
    div_ct: dict[str, int] = defaultdict(int)
    for c in enriched_paired:
        cls = c.get("class", "?")
        class_ct[cls] += 1
        if cls == "A":
            div_ct["divergent" if c.get("divergent") else "same_strand"] += 1
    for k, v in sorted(class_ct.items()):
        print(f"    {k}: {v}", file=sys.stderr)
    print(f"    class A orientation: {dict(div_ct)}", file=sys.stderr)

    print("\n=== class A IGS length distribution ===", file=sys.stderr)
    lens = [c.get("total_igs_bp", 0) for c in enriched_paired if c.get("class") == "A"]
    if lens:
        lens_sorted = sorted(lens)
        n = len(lens_sorted)
        print(f"    n={n}  min={min(lens)}  p25={lens_sorted[n // 4]}  "
              f"median={lens_sorted[n // 2]}  p75={lens_sorted[3 * n // 4]}  "
              f"max={max(lens)}", file=sys.stderr)
        print(f"    (Carey MG1655 torT-torS intergenic ~110 bp)", file=sys.stderr)

    print("\n=== unique taxa in class A ===", file=sys.stderr)
    taxa = sorted({c.get("taxon", "") for c in enriched_paired if c.get("class") == "A"})
    print(f"    {len(taxa)} distinct taxon IDs", file=sys.stderr)
    print(f"    first 10: {taxa[:10]}", file=sys.stderr)

    print(f"\nwrote {RESULTS / 'candidate_loci_enriched.tsv'}", file=sys.stderr)
    print(f"wrote {RESULTS / 'half_pair_loci_enriched.tsv'}", file=sys.stderr)
    print(f"wrote {RESULTS / 'class_A_igs.fasta'} ({len(fasta_A)} records)",
          file=sys.stderr)
    print(f"wrote {RESULTS / 'class_B_igs.fasta'} ({len(fasta_B)} records)",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
