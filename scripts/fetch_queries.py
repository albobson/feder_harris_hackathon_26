"""Fetch protein FASTAs for the pilot mmseqs search.

Writes:
    queries/torT_torS.faa        MG1655 TorT (P36664) + TorS (P39453)
    queries/phage_integrases.faa Integrase protein from each phage in Supp 1

Run with the omg_search micromamba env (PYTHONNOUSERSITE=1 recommended
while the stale user-site pandas exists):

    PYTHONNOUSERSITE=1 ~/micromamba/envs/omg_search/bin/python \
        scripts/fetch_queries.py

Set NCBI_EMAIL in your environment (Entrez asks for a contact email).
Everything is fetched from NCBI Entrez — UniProt's REST cache was
returning empty responses on this cluster (2026-07-30).
"""

from __future__ import annotations

import os
import sys
import time
from pathlib import Path

from Bio import Entrez, SeqIO

REPO = Path(__file__).resolve().parent.parent
QUERIES = REPO / "queries"

Entrez.email = os.environ.get("NCBI_EMAIL", "dustinbmullaney@gmail.com")

# ---------------------------------------------------------------------------
# 1. TorT + TorS from E. coli MG1655 (via NCBI protein)
# ---------------------------------------------------------------------------

NCBI_TORT_TORS = {
    "NP_415514.1": "torT_MG1655_periplasmic_TMAO_binding",
    "NP_415513.2": "torS_MG1655_hybrid_sensor_kinase",
}

# ---------------------------------------------------------------------------
# 2. Phage integrases from Supp 1
# ---------------------------------------------------------------------------
# One RefSeq genome per phage in Supp 1. We pull the whole phage GenBank
# record and extract the CDS whose /product contains "integrase".
# HK022 is added because it is the anchor phage in the Carey paper even
# though it is not in the Supp 1 strain table.

PHAGE_ACCESSIONS = {
    "HK022": "NC_002166",
    "lambda": "NC_001416",
    "Sf6": "NC_005344",
    "Sf101": "NC_027398",
    "BP-4795": "NC_004813",
    "CDT-1phi": "NC_009514",
    "DE3": "NC_042057",
    "GF-2": "NC_026611",
    "HK620": "NC_002730",
    "SEN34": "NC_028699",
    "YYZ-2008": "NC_011356",
}


def fetch_ncbi_protein(acc: str) -> str:
    """Fetch a single NCBI protein record as FASTA text."""
    h = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
    txt = h.read()
    h.close()
    return txt


def fetch_phage_integrase(name: str, acc: str) -> list[tuple[str, str]]:
    """Fetch a phage GenBank record and return every CDS whose product
    field mentions "integrase" or "int"."""
    handle = Entrez.efetch(db="nuccore", id=acc, rettype="gb", retmode="text")
    rec = SeqIO.read(handle, "genbank")
    handle.close()

    hits = []
    for feat in rec.features:
        if feat.type != "CDS":
            continue
        product = " ".join(feat.qualifiers.get("product", [""])).lower()
        gene = " ".join(feat.qualifiers.get("gene", [""])).lower()
        if "integrase" in product or gene in {"int", "inta"}:
            aa = feat.qualifiers.get("translation", [None])[0]
            if aa is None:
                continue
            locus = feat.qualifiers.get("locus_tag", [feat.qualifiers.get("gene", ["cds"])[0]])[0]
            hits.append((f"{name}|{acc}|{locus}|{product or gene}", aa))
    return hits


def main() -> int:
    QUERIES.mkdir(parents=True, exist_ok=True)

    tor_out = QUERIES / "torT_torS.faa"
    print(f"[torT/torS] writing {tor_out}", file=sys.stderr)
    with tor_out.open("w") as fh:
        for acc, label in NCBI_TORT_TORS.items():
            fasta = fetch_ncbi_protein(acc)
            lines = fasta.splitlines()
            if not lines or not lines[0].startswith(">"):
                print(f"    !! empty response for {acc}", file=sys.stderr)
                continue
            seq_lines = [ln for ln in lines[1:] if ln]
            fh.write(f">{acc}|{label}\n")
            fh.write("\n".join(seq_lines).rstrip() + "\n")
            print(f"    {acc} {label} ({sum(len(s) for s in seq_lines)} aa)", file=sys.stderr)
            time.sleep(0.4)

    int_out = QUERIES / "phage_integrases.faa"
    print(f"[integrases] writing {int_out}", file=sys.stderr)
    with int_out.open("w") as fh:
        for name, acc in PHAGE_ACCESSIONS.items():
            hits = fetch_phage_integrase(name, acc)
            if not hits:
                print(f"    !! no integrase annotation found for {name} ({acc})", file=sys.stderr)
                continue
            for header, aa in hits:
                fh.write(f">{header}\n{aa}\n")
                print(f"    {header} ({len(aa)} aa)", file=sys.stderr)
            time.sleep(0.4)

    return 0


if __name__ == "__main__":
    sys.exit(main())
