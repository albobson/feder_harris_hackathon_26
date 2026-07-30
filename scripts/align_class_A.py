"""Step 5 finalization: MSA the class-A divergent-promoter IGS with the
E. coli MG1655 torT-torS intergenic as a reference anchor.

Input:  results/class_A_igs.fasta            (from enrich_candidates.py)
Output: results/class_A_igs_filtered.fasta   (with MG1655 reference,
                                              length-filtered)
        results/class_A_igs_aligned.fasta    (MAFFT MSA)

Run:
    ml MAFFT/7.526-GCC-13.2.0-with-extensions
    PYTHONNOUSERSITE=1 ~/micromamba/envs/omg_search/bin/python \\
        scripts/align_class_A.py
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from Bio import Entrez, SeqIO

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
Entrez.email = os.environ.get("NCBI_EMAIL", "dustinbmullaney@gmail.com")

MIN_LEN = 50
MAX_LEN = 200


def fetch_gene_coords(gene: str) -> tuple[str, int, int, str]:
    """Return (chrom_accession, start, end, strand) for a gene in MG1655."""
    h = Entrez.esearch(
        db="gene",
        term=f'{gene}[Gene Name] AND "Escherichia coli str. K-12 substr. MG1655"[Organism]',
        retmax=3,
    )
    r = Entrez.read(h)
    h.close()
    if not r["IdList"]:
        raise RuntimeError(f"no gene id for {gene}")
    h = Entrez.esummary(db="gene", id=r["IdList"][0])
    s = Entrez.read(h)
    h.close()
    doc = s["DocumentSummarySet"]["DocumentSummary"][0]
    loc = doc["GenomicInfo"][0]
    a = int(loc["ChrStart"]) + 1  # NCBI esummary uses 0-based
    b = int(loc["ChrStop"]) + 1
    strand = "+" if a < b else "-"
    return loc["ChrAccVer"], min(a, b), max(a, b), strand


def fetch_mg1655_torT_torS_igs() -> tuple[str, str]:
    """Fetch the torT-torS intergenic region from E. coli K-12 MG1655
    (RefSeq NC_000913). Divergent: torS on minus (1053434..1056178),
    torT on plus (1056261..1057289), intergenic is 82 bp (1056179..1056260).
    Falls back to hardcoded coords if NCBI Gene lookup hiccups."""
    try:
        _, torT_a, torT_b, _ = fetch_gene_coords("torT")
        chrom, torS_a, torS_b, _ = fetch_gene_coords("torS")
    except Exception as e:
        print(f"    gene lookup failed ({e}); using hardcoded MG1655 coords",
              file=sys.stderr)
        chrom = "NC_000913.3"
        torT_a, torT_b = 1056261, 1057289  # torT plus strand
        torS_a, torS_b = 1053434, 1056178  # torS minus strand

    # Intergenic sits between the closer ends
    if torT_a > torS_b:
        ig_start, ig_end = torS_b + 1, torT_a - 1
    elif torS_a > torT_b:
        ig_start, ig_end = torT_b + 1, torS_a - 1
    else:
        raise RuntimeError(
            f"torT ({torT_a}..{torT_b}) and torS ({torS_a}..{torS_b}) overlap"
        )

    # Retry the fetch a couple times — NCBI has transient failures
    seq = ""
    last_err: Exception | None = None
    for attempt in range(3):
        try:
            h = Entrez.efetch(
                db="nuccore", id=chrom, rettype="fasta", retmode="text",
                seq_start=ig_start, seq_stop=ig_end, strand=1,
            )
            fasta = h.read().splitlines()
            h.close()
            seq = "".join(fasta[1:])
            if seq:
                break
        except Exception as e:
            last_err = e
    if not seq:
        raise RuntimeError(f"could not fetch MG1655 intergenic: {last_err}")

    header = (f"MG1655_torT_torS_intergenic|{chrom}|"
              f"{ig_start}:{ig_end} len={len(seq)}")
    return header, seq


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records = []
    header = None
    parts: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(parts)))
            header = line[1:]
            parts = []
        else:
            parts.append(line.strip())
    if header is not None:
        records.append((header, "".join(parts)))
    return records


def write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    with path.open("w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def consensus(aligned: list[tuple[str, str]]) -> str:
    """Simple majority-rule consensus over an alignment."""
    if not aligned:
        return ""
    L = len(aligned[0][1])
    cons = []
    for i in range(L):
        col = [s[i] for _, s in aligned]
        counts: dict[str, int] = {}
        for c in col:
            if c == "-":
                continue
            counts[c] = counts.get(c, 0) + 1
        if not counts:
            cons.append("-")
            continue
        best = max(counts.items(), key=lambda kv: kv[1])
        freq = best[1] / len(aligned)
        cons.append(best[0].upper() if freq >= 0.5 else best[0].lower())
    return "".join(cons)


def main() -> int:
    src = RESULTS / "class_A_igs.fasta"
    if not src.exists():
        print(f"missing input: {src}", file=sys.stderr)
        return 1

    records = read_fasta(src)
    print(f"loaded {len(records)} class-A IGS records", file=sys.stderr)

    # Fetch MG1655 reference and prepend
    print("fetching MG1655 torT-torS intergenic from NCBI...", file=sys.stderr)
    try:
        ref_header, ref_seq = fetch_mg1655_torT_torS_igs()
        print(f"    got MG1655 reference: {len(ref_seq)} bp", file=sys.stderr)
    except Exception as e:
        print(f"    !! failed to fetch MG1655: {e}; proceeding without",
              file=sys.stderr)
        ref_header, ref_seq = None, None

    filtered = [(h, s) for h, s in records if MIN_LEN <= len(s) <= MAX_LEN]
    print(f"length filter {MIN_LEN}-{MAX_LEN} bp: kept {len(filtered)}/{len(records)}",
          file=sys.stderr)

    if ref_header:
        filtered = [(ref_header, ref_seq)] + filtered

    filt_out = RESULTS / "class_A_igs_filtered.fasta"
    write_fasta(filtered, filt_out)
    print(f"wrote {filt_out}", file=sys.stderr)

    aln_out = RESULTS / "class_A_igs_aligned.fasta"
    print("running mafft --auto...", file=sys.stderr)
    result = subprocess.run(
        ["mafft", "--auto", "--reorder", str(filt_out)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print("mafft failed:", file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        return 1
    aln_out.write_text(result.stdout)
    print(f"wrote {aln_out}", file=sys.stderr)

    aligned = read_fasta(aln_out)
    cons = consensus(aligned)
    L = len(cons)
    print(f"\naligned {len(aligned)} seqs, alignment length = {L}",
          file=sys.stderr)
    coverage = sum(1 for c in cons if c != "-") / L
    conserved = sum(1 for c in cons if c.isupper()) / L
    print(f"    non-gap consensus coverage: {coverage:.1%}",
          file=sys.stderr)
    print(f"    columns with >=50% majority residue: {conserved:.1%}",
          file=sys.stderr)

    print("\n=== consensus (upper = >=50% majority) ===", file=sys.stderr)
    for i in range(0, L, 80):
        print(f"  {i+1:>4}  {cons[i:i + 80]}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
