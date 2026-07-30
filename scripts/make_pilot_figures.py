"""Generate figures for the pilot results summary.

Writes:
    figures/fig1_igs_length_distribution.png
    figures/fig2_class_breakdown.png
    figures/fig3_msa_conservation.png

Run:
    PYTHONNOUSERSITE=1 ~/micromamba/envs/omg_search/bin/python \\
        scripts/make_pilot_figures.py
"""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
FIGURES = REPO / "figures"
FIGURES.mkdir(exist_ok=True)

MG1655_IGS_BP = 82  # actual torT-torS intergenic in E. coli K-12 MG1655


def read_tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records = []
    h, s = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if h is not None:
                records.append((h, "".join(s)))
            h, s = line[1:], []
        else:
            s.append(line.strip())
    if h is not None:
        records.append((h, "".join(s)))
    return records


def fig_igs_length_distribution(cands: list[dict]) -> None:
    lens = [int(c["total_igs_bp"]) for c in cands
            if c["class"] == "A" and c["total_igs_bp"]]
    fig, ax = plt.subplots(figsize=(6.4, 4.0), dpi=140)
    ax.hist(lens, bins=range(0, 260, 10), color="#4a6fa5",
            edgecolor="white", alpha=0.9)
    ax.axvline(MG1655_IGS_BP, color="#d7263d", linestyle="--", linewidth=1.5,
               label=f"E. coli MG1655 (82 bp)")
    median = sorted(lens)[len(lens) // 2]
    ax.axvline(median, color="#404040", linestyle=":", linewidth=1.2,
               label=f"OMG median ({median} bp)")
    ax.set_xlabel("torT/torS intergenic length (bp)")
    ax.set_ylabel("class-A candidate loci")
    ax.set_title(f"OMG class-A intergenic lengths match E. coli reference\n"
                 f"(n={len(lens)} contigs, 52 distinct taxa)")
    ax.legend(frameon=False, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out = FIGURES / "fig1_igs_length_distribution.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")


def fig_class_breakdown(cands: list[dict], hp: list[dict]) -> None:
    class_ct = Counter(c["class"] for c in cands)
    divergent_ct = Counter()
    for c in cands:
        if c["class"] == "A":
            divergent_ct["divergent" if c["divergent"] == "True" else "same_strand"] += 1

    hp_by_fam = Counter(r["integrase_family"] for r in hp)

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0), dpi=140)

    # --- Left panel: class breakdown of paired candidates ---
    order = ["A", "C_small_gap", "C"]
    labels = ["A: adjacent\n(baseline)", "C: gap 2-15\n(no phage)", "C: gap >15\n(rearranged?)"]
    counts = [class_ct.get(k, 0) for k in order]
    colors = ["#3f8f4a", "#c9a03a", "#c26b2b"]
    axes[0].bar(labels, counts, color=colors, edgecolor="white")
    for i, v in enumerate(counts):
        axes[0].text(i, v + 0.6, str(v), ha="center", fontsize=11)
    axes[0].set_title("Paired candidates (torT + torS on\nsame contig, n=56)")
    axes[0].set_ylabel("contigs")
    axes[0].spines[["top", "right"]].set_visible(False)
    axes[0].set_ylim(0, max(counts) * 1.15)

    # --- Middle panel: orientation of class A ---
    ori_labels = ["divergent\n(Carey-like)", "same strand"]
    ori_counts = [divergent_ct.get("divergent", 0),
                  divergent_ct.get("same_strand", 0)]
    axes[1].bar(ori_labels, ori_counts, color=["#3f8f4a", "#a8a8a8"],
                edgecolor="white")
    for i, v in enumerate(ori_counts):
        axes[1].text(i, v + 0.8, str(v), ha="center", fontsize=11)
    axes[1].set_title("Class A orientation\n(50/52 match Carey architecture)")
    axes[1].set_ylabel("contigs")
    axes[1].spines[["top", "right"]].set_visible(False)
    axes[1].set_ylim(0, max(ori_counts) * 1.15)

    # --- Right panel: half-pair (B') by integrase family ---
    fams = sorted(hp_by_fam, key=lambda k: -hp_by_fam[k])
    counts_hp = [hp_by_fam[f] for f in fams]
    bars = axes[2].barh(fams[::-1], counts_hp[::-1], color="#4a6fa5",
                        edgecolor="white")
    for bar, v in zip(bars, counts_hp[::-1]):
        axes[2].text(v + 0.15, bar.get_y() + bar.get_height() / 2, str(v),
                     va="center", fontsize=10)
    axes[2].set_title("Class B': torS + adjacent phage integrase\n"
                      "(all n=51 records are torS-anchored)")
    axes[2].set_xlabel("half-pair records")
    axes[2].spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    out = FIGURES / "fig2_class_breakdown.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")


def fig_msa_conservation(aligned: list[tuple[str, str]]) -> None:
    """Per-column conservation (frequency of most common non-gap base)
    across the class-A MSA, with a horizontal band marking where the
    MG1655 reference aligns."""
    seqs = [s for _, s in aligned]
    L = len(seqs[0])
    freqs = []
    coverage = []
    for i in range(L):
        col = [row[i] for row in seqs]
        non_gap = [c for c in col if c != "-"]
        if not non_gap:
            freqs.append(0.0)
            coverage.append(0)
            continue
        counts = Counter(non_gap)
        freqs.append(counts.most_common(1)[0][1] / len(seqs))
        coverage.append(len(non_gap))

    mg = next((s for h, s in aligned if "MG1655" in h), None)
    mg_lo = len(mg) - len(mg.lstrip("-")) if mg else 0
    mg_hi = len(mg) - len(mg.rstrip("-")) if mg else 0
    mg_hi = len(mg) - mg_hi if mg else L

    fig, ax = plt.subplots(figsize=(9.0, 3.6), dpi=140)
    ax.fill_between(range(L), 0, freqs, color="#4a6fa5", alpha=0.8,
                    linewidth=0)
    if mg:
        ax.axvspan(mg_lo, mg_hi, color="#d7263d", alpha=0.12,
                   label=f"MG1655 aligned span (cols {mg_lo + 1}-{mg_hi})")
    ax.set_xlabel("MSA column")
    ax.set_ylabel("majority-base frequency")
    ax.set_title(f"Column-wise conservation of class-A intergenic MSA "
                 f"(n={len(seqs)} sequences, {L} cols)")
    ax.set_ylim(0, 1.05)
    ax.set_xlim(0, L)
    if mg:
        ax.legend(loc="upper right", frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out = FIGURES / "fig3_msa_conservation.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")


def main() -> int:
    cands = read_tsv(RESULTS / "candidate_loci_enriched.tsv")
    hp = read_tsv(RESULTS / "half_pair_loci_enriched.tsv")
    aligned = read_fasta(RESULTS / "class_A_igs_aligned.fasta")

    fig_igs_length_distribution(cands)
    fig_class_breakdown(cands, hp)
    fig_msa_conservation(aligned)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
