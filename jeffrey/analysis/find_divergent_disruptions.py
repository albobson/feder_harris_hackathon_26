import os, re, json

# Repo-relative paths.
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(os.path.join(HERE, "..", "data"))
GDIR = os.path.join(DATA_DIR, "genomes")

# For each species, one "baseline" strain used to define the set of divergent gene
# pairs with a short (shared-regulatory-region-sized) intergenic spacer. Other
# strains of the same species are then checked for the same gene pair having been
# pulled apart by an intervening insertion. Picked as the well-known type/lab strain
# for each species so the baseline itself is unlikely to carry an unusual insertion.
REFERENCE_STRAIN = {
    "Ecoli": "GCF_000005845.2",         # K-12 MG1655
    "Saureus": "GCF_000013425.1",       # NCTC 8325
    "Paeruginosa": "GCF_000006765.1",   # PAO1
    "Vcholerae": "GCF_055389885.1",     # N16961-like / RIMD 2214379
    "Mtuberculosis": "GCF_000195955.2", # H37Rv
    "Bsubtilis": "GCF_000009045.1",     # 168
    "Spyogenes": "GCF_002055535.1",     # NCTC8198
    "Bfragilis": "GCF_000025985.1",     # NCTC 9343
    "Cjejuni": "GCF_000009085.1",       # NCTC 11168
}

ATTR_RE = re.compile(r'([^=;]+)=([^;]*)')

def parse_attrs(s):
    # GFF3 column 9 is "key=value;key=value;...".
    return dict(ATTR_RE.findall(s))

def load_gff(path):
    """Parse a GFF3 file into a flat list of gene records with strand/coords/symbol,
    and attach each gene's CDS product description (used later to recognize
    integrase/transposase/phage genes sitting inside an expanded intergenic region).
    Also returns each replicon's full length (from the ##sequence-region pragma),
    needed to compute distances correctly across the origin of a circular chromosome."""
    genes = []
    products = {}  # gene ID -> product string, collected from CDS features
    seq_lengths = {}  # seqid -> replicon length, e.g. "##sequence-region NZ_CP... 1 4459449"
    with open(path) as f:
        for line in f:
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4:
                    seq_lengths[parts[1]] = int(parts[3])
                continue
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, frame, attrs = cols
            a = parse_attrs(attrs)
            if ftype == "gene":
                # RefSeq/PGAP annotation puts the gene symbol in "gene=" (fallback "Name=").
                # This symbol is what lets us match the "same" gene across different
                # strains' GFFs, since locus_tags differ per assembly.
                symbol = a.get("gene") or a.get("Name")
                locus_tag = a.get("locus_tag")
                gid = a.get("ID")
                genes.append({
                    "seqid": seqid, "start": int(start), "end": int(end),
                    "strand": strand, "symbol": symbol, "locus_tag": locus_tag,
                    "id": gid, "product": None,
                })
            elif ftype == "CDS":
                # In NCBI prokaryotic GFF3, CDS Parent= points directly at the gene ID.
                parent = a.get("Parent")
                prod = a.get("product")
                if parent and prod:
                    products[parent] = prod
    for g in genes:
        if g["id"] in products:
            g["product"] = products[g["id"]]
    return genes, seq_lengths

def divergent_pairs(genes, max_len=500):
    """Find pairs of immediately-adjacent genes (no gene between them) transcribed
    divergently -- left gene on '-' strand, right gene on '+' strand -- so both
    promoters face into, and potentially share regulatory features within, the
    intergenic spacer between them. Restricted to spacers <= max_len bp, the size
    range typical of a real shared promoter/operator region (c.f. torS/torT ~130bp).
    Keyed by the pair of gene symbols so it can be looked up in other genomes."""
    by_seq = {}
    for g in genes:
        by_seq.setdefault(g["seqid"], []).append(g)
    pairs = {}
    for seqid, glist in by_seq.items():
        glist.sort(key=lambda g: g["start"])
        for i in range(len(glist) - 1):
            a, b = glist[i], glist[i + 1]
            if a["strand"] == "-" and b["strand"] == "+":
                length = b["start"] - a["end"] - 1
                if a["symbol"] and b["symbol"] and 0 <= length <= max_len:
                    key = tuple(sorted([a["symbol"], b["symbol"]]))
                    pairs[key] = {"seqid": seqid, "length": length,
                                  "a": a["symbol"], "b": b["symbol"]}
    return pairs

def gene_lookup(genes):
    # symbol -> list of gene records (usually length 1; a list in case a symbol repeats).
    d = {}
    for g in genes:
        if g["symbol"]:
            d.setdefault(g["symbol"], []).append(g)
    return d

def find_span(lookup, symA, symB, seq_lengths):
    """In one genome, locate symA and symB and, if they are still divergently
    oriented on the same replicon, measure the distance between them going the
    SHORT way around (chromosomes are circular, so we compute the wraparound
    distance too -- otherwise a pair that sits near wherever this particular
    assembly's linear coordinates happen to start/end looks artifactually
    "separated by the whole genome" even though nothing was inserted). In a
    strain with no insertion this should match the reference spacer; in a strain
    where a mobile element has landed between them, the true (possibly wrapped)
    distance will be much larger."""
    if symA not in lookup or symB not in lookup:
        return None
    best = None
    for gA in lookup[symA]:
        for gB in lookup[symB]:
            if gA["seqid"] != gB["seqid"]:
                continue
            L = seq_lengths.get(gA["seqid"])
            # Try both assignments of which gene is "left" (the divergent pair's
            # orientation, not raw coordinate order, decides which is which).
            for left, right in ((gA, gB), (gB, gA)):
                if left["strand"] != "-" or right["strand"] != "+":
                    continue
                if L:
                    span = (right["start"] - left["end"] - 1) % L
                else:
                    span = right["start"] - left["end"] - 1
                cand = {"seqid": gA["seqid"], "span": span, "left": left,
                         "right": right, "length": L,
                         "wrapped": right["start"] <= left["end"]}
                if best is None or cand["span"] < best["span"]:
                    best = cand
    return best

def intervening_products(genes, seqid, left_end, right_start, length, wrapped):
    # Genes fully contained within the (now-expanded) intergenic window -- i.e. the
    # candidate inserted element's own genes -- so we can eyeball what they encode.
    # If the true short path wraps past the origin, the window is the union of
    # (left_end, replicon end] and [replicon start, right_start).
    prods = []
    for g in genes:
        if g["seqid"] != seqid:
            continue
        if wrapped:
            inside = g["start"] > left_end or g["end"] < right_start
        else:
            inside = g["start"] >= left_end and g["end"] <= right_start
        if inside:
            prods.append((g["symbol"] or g["locus_tag"], g["product"]))
    return prods

def run_full_scan():
    """The original entry point: for every downloaded species, find all short-spacer
    divergent pairs in the reference strain and check every other strain for
    disruption. Wrapped in a function (rather than bare module-level code) so other
    scripts can `import find_divergent_disruptions as fdd` and reuse its parsing
    helpers (load_gff, divergent_pairs, gene_lookup, find_span, intervening_products)
    without re-running this whole scan as a side effect of the import."""
    results = {}

    for species_dir in sorted(os.listdir(GDIR)):
        sp_path = os.path.join(GDIR, species_dir)
        if not os.path.isdir(sp_path):
            continue
        ref_acc = REFERENCE_STRAIN.get(species_dir)
        accs = sorted(os.listdir(sp_path))
        if ref_acc not in accs:
            ref_acc = accs[0]  # fallback if the intended reference wasn't downloaded

        genome_genes = {}
        genome_seqlens = {}
        for acc in accs:
            gff = os.path.join(sp_path, acc, "genomic.gff")
            if os.path.exists(gff):
                genes, seq_lengths = load_gff(gff)
                genome_genes[acc] = genes
                genome_seqlens[acc] = seq_lengths

        ref_genes = genome_genes.get(ref_acc)
        if not ref_genes:
            continue
        ref_pairs = divergent_pairs(ref_genes, max_len=500)
        print(f"{species_dir}: reference {ref_acc}, {len(ref_pairs)} divergent pairs (<=500bp spacer) found")

        species_hits = []
        for acc, genes in genome_genes.items():
            if acc == ref_acc:
                continue
            lookup = gene_lookup(genes)
            seq_lengths = genome_seqlens[acc]
            for (symA, symB), refinfo in ref_pairs.items():
                span_info = find_span(lookup, symA, symB, seq_lengths)
                if span_info is None:
                    continue
                # Flag as a candidate disruption if the spacer blew up far beyond the
                # reference size (>5x, and an absolute floor of 800bp so noise in small
                # reference spacers doesn't trigger false positives).
                if span_info["span"] > 5 * max(refinfo["length"], 50) and span_info["span"] > 800:
                    prods = intervening_products(genes, span_info["seqid"],
                                                  span_info["left"]["end"], span_info["right"]["start"],
                                                  span_info["length"], span_info["wrapped"])
                    species_hits.append({
                        "genome": acc, "pair": (symA, symB),
                        "ref_length": refinfo["length"], "expanded_length": span_info["span"],
                        "n_intervening_genes": len(prods),
                        "intervening_products": prods,
                    })
        results[species_dir] = {"reference": ref_acc, "n_ref_pairs": len(ref_pairs), "hits": species_hits}

    with open(os.path.join(DATA_DIR, "divergent_disruption_results.json"), "w") as f:
        json.dump(results, f, indent=2)

    print("\n=== SUMMARY ===")
    for sp, r in results.items():
        print(f"{sp}: {len(r['hits'])} candidate disruption(s) across {r['n_ref_pairs']} reference divergent pairs")
    return results


if __name__ == "__main__":
    run_full_scan()
