"""Sequence-based ortholog mapping from our reference genomes to STRING protein
IDs, replacing the brittle gene-symbol join.

WHY: STRING's alias tables carry gene symbols for only some organisms. For
V. cholerae, S. pyogenes and B. fragilis the symbol join mapped just 1-6% of
divergent pairs, because STRING's reference strain is annotated with locus tags
(and, for V. cholerae, 2000-era GenBank protein accessions like AAF93179.1)
that share no namespace with our RefSeq WP_ accessions. Since these are
within-species comparisons, protein sequences are near-identical and homology
search maps them reliably.

Only proteins belonging to genes that appear in a divergent pair are searched,
which keeps the BLAST work small.

Mapping is one-directional best-hit (by bitscore) under strict thresholds.
That is defensible here because query and subject are the SAME SPECIES, so the
true ortholog is typically a >95%-identity, unambiguous top hit; reciprocal
best hits would be the stricter choice for cross-species work.
"""
import glob, gzip, json, os, re, shutil, subprocess, sys
import find_divergent_disruptions as fdd

DATA_DIR = fdd.DATA_DIR
GDIR = fdd.GDIR
SDIR = os.path.join(DATA_DIR, "string")
PROTDIR = os.path.join(SDIR, "proteomes")
REFPROT = os.path.join(DATA_DIR, "refprot")
WORK = os.path.join(DATA_DIR, "blast_work")

MIN_PIDENT = 70.0   # same-species: expect very high identity
MIN_QCOV = 0.70     # guard against short partial alignments
EVALUE = "1e-20"
# Respect the slot count SGE granted (NSLOTS) rather than hardcoding threads,
# so this stays polite on a shared cluster. Falls back to 1 when run outside SGE.
THREADS = os.environ.get("NSLOTS", "1")

ATTR_RE = re.compile(r'([^=;]+)=([^;]*)')


def gene_to_protein(gff_path):
    """Map gene symbol -> protein_id, via CDS Parent= pointing at the gene ID."""
    gid_to_symbol, sym_to_prot = {}, {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            ftype, attrs = cols[2], cols[8]
            a = dict(ATTR_RE.findall(attrs))
            if ftype == "gene":
                sym = a.get("gene") or a.get("Name")
                if a.get("ID") and sym:
                    gid_to_symbol[a["ID"]] = sym
            elif ftype == "CDS":
                parent, pid = a.get("Parent"), a.get("protein_id")
                if parent in gid_to_symbol and pid:
                    sym_to_prot.setdefault(gid_to_symbol[parent], pid)
    return sym_to_prot


def read_fasta(path):
    opener = gzip.open if path.endswith(".gz") else open
    seqs, name, buf = {}, None, []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def main():
    os.makedirs(WORK, exist_ok=True)
    taxids = {}
    with open(os.path.join(SDIR, "taxids.tsv")) as f:
        for line in f:
            name, taxid = line.strip().split("\t")
            taxids[name] = taxid

    ortholog_map = {}
    stats = {}

    for species in sorted(taxids):
        taxid = taxids[species]
        sp_path = os.path.join(GDIR, species)
        ref_acc = fdd.REFERENCE_STRAIN.get(species)
        accs = sorted(os.listdir(sp_path))
        if ref_acc not in accs:
            ref_acc = accs[0]
        gff = os.path.join(sp_path, ref_acc, "genomic.gff")

        ref_genes, _ = fdd.load_gff(gff)
        pairs = fdd.divergent_pairs(ref_genes, max_len=500)
        needed_syms = set()
        for a, b in pairs:
            needed_syms.add(a)
            needed_syms.add(b)

        sym_to_prot = gene_to_protein(gff)
        # proteins we actually need to search
        prot_to_sym = {}
        for sym in needed_syms:
            pid = sym_to_prot.get(sym)
            if pid:
                prot_to_sym.setdefault(pid, sym)

        # --- write query subset FASTA
        ref_faa = os.path.join(REFPROT, f"{species}.faa")
        all_seqs = read_fasta(ref_faa)
        qpath = os.path.join(WORK, f"{species}.query.faa")
        n_q = 0
        with open(qpath, "w") as out:
            for pid, sym in prot_to_sym.items():
                seq = all_seqs.get(pid)
                if not seq:
                    continue
                out.write(f">{pid}\n{seq}\n")
                n_q += 1

        # --- build STRING blast db (gunzip first; makeblastdb can't read .gz)
        dbfa = os.path.join(WORK, f"{taxid}.string.faa")
        if not os.path.exists(dbfa):
            with gzip.open(os.path.join(PROTDIR, f"{taxid}.fa.gz"), "rt") as fi, open(dbfa, "w") as fo:
                shutil.copyfileobj(fi, fo)
        if not os.path.exists(dbfa + ".phr"):
            subprocess.run(["makeblastdb", "-in", dbfa, "-dbtype", "prot", "-out", dbfa],
                           check=True, stdout=subprocess.DEVNULL)

        # --- blastp
        tsv = os.path.join(WORK, f"{species}.hits.tsv")
        if not os.path.exists(tsv):
            subprocess.run([
                "blastp", "-query", qpath, "-db", dbfa, "-out", tsv,
                "-outfmt", "6 qseqid sseqid pident length qlen bitscore evalue",
                "-evalue", EVALUE, "-max_target_seqs", "5", "-num_threads", THREADS,
            ], check=True)

        # --- best hit per query, under thresholds
        best = {}
        with open(tsv) as f:
            for line in f:
                q, s, pident, length, qlen, bitscore, ev = line.rstrip("\n").split("\t")
                pident, length, qlen, bitscore = float(pident), int(length), int(qlen), float(bitscore)
                if pident < MIN_PIDENT or (length / qlen) < MIN_QCOV:
                    continue
                if q not in best or bitscore > best[q][1]:
                    best[q] = (s, bitscore)

        sym_map = {}
        for pid, (sid, _) in best.items():
            sym = prot_to_sym.get(pid)
            if sym:
                sym_map[sym] = sid
        ortholog_map[species] = sym_map

        stats[species] = {
            "taxid": taxid, "genes_in_pairs": len(needed_syms),
            "queried": n_q, "mapped": len(sym_map),
            "rate": round(len(sym_map) / len(needed_syms), 3) if needed_syms else None,
        }
        print(f"{species:<14} genes_in_pairs={len(needed_syms):<5} queried={n_q:<5} "
              f"ortholog-mapped={len(sym_map):<5} ({stats[species]['rate']:.0%})")

    with open(os.path.join(DATA_DIR, "ortholog_map.json"), "w") as f:
        json.dump({"stats": stats, "map": ortholog_map,
                   "thresholds": {"min_pident": MIN_PIDENT, "min_qcov": MIN_QCOV,
                                  "evalue": EVALUE, "method": "one-directional best hit by bitscore"}},
                  f, indent=2)
    print("\nwrote ortholog_map.json")


if __name__ == "__main__":
    main()
