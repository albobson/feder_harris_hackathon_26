"""Does the STRING coupling score reproduce our hand-verified calls?

Parts 2 and 4 produced a labelled set by literature checking: pairs confirmed
functionally coupled, and pairs explicitly confirmed NOT coupled (regulator
present but its real target is elsewhere, or functions unrelated). That is a
calibration set. If the automated score separates them, it can be trusted at
scale; if not, the annotation layer is not yet fixed.

Also quantifies the circularity trap: the `neighborhood` channel should be
high for BOTH labels (all these pairs are adjacent), demonstrating that it
carries adjacency rather than function and had to be excluded.
"""
import csv, os
import find_divergent_disruptions as fdd

DATA_DIR = fdd.DATA_DIR

# (species, geneA, geneB, label, note). label: 1 = verified coupled, 0 = verified NOT coupled
LABELS = [
    # --- verified functionally coupled (Part 4 curated + confirmed heuristic hits)
    ("Ecoli", "torS", "torT", 1, "same TMAO two-component pathway (the paper's case)"),
    ("Ecoli", "metE", "metR", 1, "MetR activates metE"),
    ("Vcholerae", "metE", "metR", 1, "MetR activates metE"),
    ("Ecoli", "araC", "araB", 1, "AraC activates araBAD"),
    ("Ecoli", "soxR", "soxS", 1, "SoxR activates soxS"),
    ("Ecoli", "acrR", "acrA", 1, "AcrR represses acrAB"),
    ("Ecoli", "cpxP", "cpxR", 1, "CpxP accessory to CpxRA"),
    ("Bsubtilis", "rocD", "rocR", 1, "RocR activates rocDEF operon"),
    ("Bsubtilis", "yyaT", "yybA", 1, "YybA regulon includes yyaT (predicted)"),
    # --- verified NOT coupled (Part 2 literature checks + Part 4 refutations)
    ("Ecoli", "fabR", "sthA", 0, "FabR regulon is fabA/fabB, not sthA"),
    ("Ecoli", "btsS", "mlrA", 0, "MlrA's real target is csgD"),
    ("Bsubtilis", "purT", "ybfI", 0, "YbfI function uncharacterised"),
    ("Ecoli", "flk", "pdxB", 0, "documented coincidental divergent promoter"),
    ("Ecoli", "pepE", "rluF", 0, "unrelated functions"),
    ("Ecoli", "dusC", "yohJ", 0, "unrelated functions"),
    ("Ecoli", "ybjP", "ybjQ", 0, "unrelated functions"),
    ("Ecoli", "dnaB", "qorA", 0, "unrelated functions"),
    ("Ecoli", "tonB", "yciI", 0, "no known link (yciI uncharacterised)"),
    ("Ecoli", "btsT", "tsr", 0, "no direct connection"),
]
AMBIGUOUS = [("Ecoli", "wrbA", "ymdF", "both RpoS regulon members but no direct link")]


def load_scores():
    scores = {}
    with open(os.path.join(DATA_DIR, "coupling_scores.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            key = (row["species"], tuple(sorted((row["geneA"], row["geneB"]))))
            scores[key] = row
    return scores


def fnum(v):
    return None if v in (None, "") else float(v)


def auc(pos, neg):
    """Rank-based AUC (Mann-Whitney): P(random positive scores above random negative)."""
    if not pos or not neg:
        return None
    wins = ties = 0
    for p in pos:
        for n in neg:
            if p > n:
                wins += 1
            elif p == n:
                ties += 1
    return (wins + 0.5 * ties) / (len(pos) * len(neg))


def median(xs):
    if not xs:
        return None
    s = sorted(xs)
    m = len(s) // 2
    return s[m] if len(s) % 2 else (s[m - 1] + s[m]) / 2


def main():
    scores = load_scores()
    print("=" * 100)
    print("VALIDATION: STRING non-positional coupling score vs hand-verified labels")
    print("=" * 100)
    print(f"{'label':<9} {'species':<13} {'pair':<16} {'coupling':>9} {'nbrhood':>8} "
          f"{'exper':>6} {'datab':>6} {'textm':>6}  note")
    print("-" * 100)

    pos, neg = [], []
    missing = []
    for species, a, b, label, note in LABELS:
        row = scores.get((species, tuple(sorted((a, b)))))
        tag = "COUPLED" if label else "not"
        if row is None:
            print(f"{tag:<9} {species:<13} {a+'/'+b:<16} {'PAIR NOT FOUND':>9}  ({note})")
            missing.append((species, a, b))
            continue
        s = fnum(row["score_primary"])
        if s is None:
            print(f"{tag:<9} {species:<13} {a+'/'+b:<16} {'unmapped':>9}  ({note})")
            missing.append((species, a, b))
            continue
        nb = fnum(row["neighborhood"]) or 0
        print(f"{tag:<9} {species:<13} {a+'/'+b:<16} {s:>9.0f} {nb:>8.0f} "
              f"{fnum(row['experimental']) or 0:>6.0f} {fnum(row['database']) or 0:>6.0f} "
              f"{fnum(row['textmining']) or 0:>6.0f}  {note}")
        (pos if label else neg).append(s)

    print()
    for species, a, b, note in AMBIGUOUS:
        row = scores.get((species, tuple(sorted((a, b)))))
        s = fnum(row["score_primary"]) if row else None
        print(f"{'ambig':<9} {species:<13} {a+'/'+b:<16} "
              f"{('%.0f' % s) if s is not None else 'unmapped':>9}  {note}")

    print()
    print("=" * 100)
    print("SEPARATION")
    print("=" * 100)
    print(f"verified COUPLED     n={len(pos):<3} median coupling score = {median(pos)}")
    print(f"verified NOT coupled n={len(neg):<3} median coupling score = {median(neg)}")
    a = auc(pos, neg)
    print(f"AUC (P[coupled scores above not-coupled]) = {a:.3f}" if a else "AUC: n/a")
    if a is not None:
        verdict = ("EXCELLENT -- score reproduces hand calls" if a >= 0.9 else
                   "GOOD -- usable with a threshold" if a >= 0.8 else
                   "MARGINAL -- better than the keyword heuristic but noisy" if a >= 0.65 else
                   "POOR -- annotation layer still not fixed")
        print(f"verdict: {verdict}")
    if missing:
        print(f"\nunscored labelled pairs ({len(missing)}): {missing}")

    # --- the circularity demonstration
    print()
    print("=" * 100)
    print("WHY `neighborhood` HAD TO BE EXCLUDED")
    print("=" * 100)
    nb_pos, nb_neg = [], []
    for species, a_, b_, label, _ in LABELS:
        row = scores.get((species, tuple(sorted((a_, b_)))))
        if not row or fnum(row["score_primary"]) is None:
            continue
        (nb_pos if label else nb_neg).append(fnum(row["neighborhood"]) or 0)
    print(f"median `neighborhood` score, verified COUPLED     = {median(nb_pos)}")
    print(f"median `neighborhood` score, verified NOT coupled = {median(nb_neg)}")
    nb_auc = auc(nb_pos, nb_neg)
    print(f"AUC using `neighborhood` alone = {nb_auc:.3f}" if nb_auc else "")
    print("If that AUC is near 0.5 while the coupling AUC is high, `neighborhood`")
    print("is tracking adjacency (which every pair here shares) rather than function")
    print("-- exactly the circularity we excluded it to avoid.")


if __name__ == "__main__":
    main()
