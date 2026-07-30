"""Step 3 of the plan: classify torT/torS/integrase hits by contig
architecture.

Input:  results/torT_torS_vs_omg10M.m8
        results/integrases_vs_omg10M.m8
Output: results/candidate_loci.tsv       (paired torT + torS candidates)
        results/half_pair_loci.tsv       (class B' — one anchor + adjacent phage)

Architecture classes (see mmseqs_search.md step 3):
    A  = torT/torS pair adjacent (|i-j| == 1) with no intervening phage.
    B  = torT/torS pair with at least one integrase hit between them.
    C  = torT/torS pair separated by many CDS with no phage signature.
    B' = one anchor + adjacent integrase, partner missing from this contig.

Notes:
- Uses cluster representatives only — a single .m8 hit represents a cluster
  that may include many source CDS. Expansion to members is a later step.
- With torS truncated at --max-seqs 5000 in the pilot search, this
  undercounts candidates on the torS side. Punting to a later rerun.

Run:
    PYTHONNOUSERSITE=1 ~/micromamba/envs/omg_search/bin/python \\
        scripts/classify_candidates.py
"""

from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"

TORT_HITS = RESULTS / "torT_torS_vs_omg10M.m8"
INT_HITS = RESULTS / "integrases_vs_omg10M.m8"

# Adjacency thresholds (CDS-index units).
# Strictly adjacent (|i-j| == 1) is the canonical Carey architecture.
# Loosen up to K to catch cases where a small insertion has occurred.
CANDIDATE_K = 15

# Column order from the mmseqs --format-output in the sbatch scripts.
M8_COLS = ["query", "target", "theader", "pident", "alnlen",
           "evalue", "bits", "qcov", "tcov"]

HEADER_RE = re.compile(r"^omg_row(\d+)_cds(\d+)$")


def parse_header(theader: str) -> tuple[int, int] | None:
    m = HEADER_RE.match(theader)
    return (int(m.group(1)), int(m.group(2))) if m else None


def classify_tort(query: str) -> str | None:
    """Return 'torT' or 'torS' from the mmseqs query header, else None."""
    q = query.lower()
    if "|tort_" in q:
        return "torT"
    if "|tors_" in q:
        return "torS"
    return None


def integrase_family(query: str) -> str:
    """Extract the phage family from the integrase query header
    (e.g. 'HK022|NC_002166|...' -> 'HK022')."""
    return query.split("|", 1)[0]


def load_hits(path: Path) -> list[dict]:
    """Load a mmseqs .m8 as a list of dicts, one per hit."""
    rows = []
    with path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        for r in reader:
            if len(r) != len(M8_COLS):
                continue
            d = dict(zip(M8_COLS, r))
            d["bits"] = float(d["bits"])
            d["evalue"] = float(d["evalue"])
            d["pident"] = float(d["pident"])
            d["qcov"] = float(d["qcov"])
            d["tcov"] = float(d["tcov"])
            key = parse_header(d["theader"])
            if key is None:
                continue
            d["row"], d["cds"] = key
            rows.append(d)
    return rows


def index_by_row(hits: list[dict], class_fn) -> dict[int, dict[str, list[dict]]]:
    """Group hits by (row, class). Returns {row: {class: [hit, hit, ...]}}."""
    out: dict[int, dict[str, list[dict]]] = defaultdict(lambda: defaultdict(list))
    for h in hits:
        cls = class_fn(h["query"])
        if cls is None:
            continue
        h["class"] = cls
        out[h["row"]][cls].append(h)
    return out


def best_hit(hits: list[dict]) -> dict:
    """Best hit by bit-score (ties broken by lower evalue)."""
    return max(hits, key=lambda h: (h["bits"], -h["evalue"]))


def main() -> int:
    tort_hits = load_hits(TORT_HITS)
    int_hits = load_hits(INT_HITS)
    print(f"loaded {len(tort_hits)} torT/torS hits, {len(int_hits)} integrase hits",
          file=sys.stderr)

    tort_by_row = index_by_row(tort_hits, classify_tort)  # {row: {'torT':[...], 'torS':[...]}}
    int_by_row = index_by_row(int_hits, integrase_family)  # {row: {family:[...]}}

    rows_torT = {r for r, v in tort_by_row.items() if "torT" in v}
    rows_torS = {r for r, v in tort_by_row.items() if "torS" in v}
    rows_int = set(int_by_row.keys())

    print(f"contig rows with torT hit: {len(rows_torT)}", file=sys.stderr)
    print(f"contig rows with torS hit: {len(rows_torS)}", file=sys.stderr)
    print(f"contig rows with any integrase hit: {len(rows_int)}", file=sys.stderr)
    both = rows_torT & rows_torS
    print(f"contig rows with BOTH torT and torS: {len(both)}", file=sys.stderr)
    print(f"    ...of which also have an integrase hit: {len(both & rows_int)}",
          file=sys.stderr)

    # ------------------------------------------------------------------
    # A/B/C — paired candidates
    # ------------------------------------------------------------------
    pair_out = RESULTS / "candidate_loci.tsv"
    pair_fields = [
        "omg_row", "class", "torT_cds", "torS_cds", "gap",
        "torT_bits", "torT_pident", "torS_bits", "torS_pident",
        "n_integrase_between", "integrase_families_between",
        "integrase_cds_between", "notes",
    ]
    class_counts: dict[str, int] = defaultdict(int)

    with pair_out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=pair_fields, delimiter="\t")
        w.writeheader()
        for row in sorted(both):
            torT_best = best_hit(tort_by_row[row]["torT"])
            torS_best = best_hit(tort_by_row[row]["torS"])
            i, j = torT_best["cds"], torS_best["cds"]
            gap = abs(i - j)
            lo, hi = min(i, j), max(i, j)

            # integrase hits that sit strictly between the two anchors
            int_between = []
            for fam, hs in int_by_row.get(row, {}).items():
                for h in hs:
                    if lo < h["cds"] < hi:
                        int_between.append((h["cds"], fam, h["bits"]))
            fams_between = sorted({f for _, f, _ in int_between})

            if gap == 0:
                cls = "?fused"
                notes = "torT and torS hits at same CDS index (possible fusion or error)"
            elif gap == 1:
                cls = "A"
                notes = ""
            elif int_between:
                cls = "B"
                notes = f"phage insertion inferred; gap={gap}"
            elif gap <= CANDIDATE_K:
                cls = "C_small_gap"
                notes = f"pair separated by {gap} CDS, no phage hit — small gap"
            else:
                cls = "C"
                notes = f"pair separated by {gap} CDS, no phage hit"

            class_counts[cls] += 1

            w.writerow({
                "omg_row": row,
                "class": cls,
                "torT_cds": i,
                "torS_cds": j,
                "gap": gap,
                "torT_bits": f'{torT_best["bits"]:.1f}',
                "torT_pident": f'{torT_best["pident"]:.2f}',
                "torS_bits": f'{torS_best["bits"]:.1f}',
                "torS_pident": f'{torS_best["pident"]:.2f}',
                "n_integrase_between": len(int_between),
                "integrase_families_between": ",".join(fams_between),
                "integrase_cds_between": ",".join(str(c) for c, _, _ in int_between),
                "notes": notes,
            })

    print(f"wrote {pair_out}", file=sys.stderr)
    print("class breakdown (paired):", dict(class_counts), file=sys.stderr)

    # ------------------------------------------------------------------
    # B' — half-pair (one anchor + adjacent integrase, partner missing)
    # ------------------------------------------------------------------
    hp_out = RESULTS / "half_pair_loci.tsv"
    hp_fields = [
        "omg_row", "anchor", "anchor_cds", "anchor_bits", "anchor_pident",
        "integrase_cds", "integrase_family", "integrase_bits", "gap_to_integrase",
    ]
    half_pair_count = 0

    with hp_out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=hp_fields, delimiter="\t")
        w.writeheader()
        for row in sorted((rows_torT ^ rows_torS) & rows_int):
            # xor: contig has exactly one of torT / torS
            anchor = "torT" if row in rows_torT else "torS"
            anchor_hit = best_hit(tort_by_row[row][anchor])
            for fam, hs in int_by_row[row].items():
                for h in hs:
                    gap = abs(anchor_hit["cds"] - h["cds"])
                    if gap == 0 or gap > CANDIDATE_K:
                        continue
                    w.writerow({
                        "omg_row": row,
                        "anchor": anchor,
                        "anchor_cds": anchor_hit["cds"],
                        "anchor_bits": f'{anchor_hit["bits"]:.1f}',
                        "anchor_pident": f'{anchor_hit["pident"]:.2f}',
                        "integrase_cds": h["cds"],
                        "integrase_family": fam,
                        "integrase_bits": f'{h["bits"]:.1f}',
                        "gap_to_integrase": gap,
                    })
                    half_pair_count += 1

    print(f"wrote {hp_out} ({half_pair_count} half-pair records)", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
