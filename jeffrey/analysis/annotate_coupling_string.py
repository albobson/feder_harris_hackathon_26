"""Annotation layer: attach a CONTINUOUS functional-coupling score to every
divergent gene pair, from STRING v12, so the selection test no longer depends
on hand-verifying a handful of pairs.

WHY STRING AND NOT KEGG: KEGG pathway co-membership misses regulator->target
relationships, which are exactly the architecture we care about. MetR is a
transcriptional activator, not a member of the methionine-biosynthesis map, so
KEGG would score metE/metR -- our strongest Part 4 hit -- as uncoupled.

THE CIRCULARITY TRAP: STRING's headline `combined_score` folds in a
`neighborhood` channel computed from CONSERVED GENOMIC ADJACENCY, plus `fusion`
and `cooccurence` which are likewise genomic-context-derived. Our pairs are
adjacent by construction, so scoring them with combined_score would reward them
for being neighbours and manufacture the very association we are trying to
test. This script therefore builds its coupling score from NON-POSITIONAL
channels only and reports the positional ones separately, so the confound can
be shown rather than assumed.

Channel usage:
  EXCLUDED always : neighborhood, fusion, cooccurence   (genomic-context)
  PRIMARY score   : experimental + database + textmining
  SENSITIVITY     : primary + coexpression. Coexpression is not positional,
                    but two genes sharing a divergent promoter are co-regulated
                    (hence co-expressed) as a consequence of the architecture
                    itself, so it is reported separately rather than trusted.
"""
import gzip, json, os, re
import find_divergent_disruptions as fdd

DATA_DIR = fdd.DATA_DIR
GDIR = fdd.GDIR
SDIR = os.path.join(DATA_DIR, "string")
STRING_VERSION = "v12.0"

MIN_SIZE, MAX_SIZE = 800, 200_000  # discrete-insertion window, as in Parts 1/3
MOBILE_KW = re.compile(
    r"integrase|transposase|phage|recombinase|capsid|tail|terminase|portal|"
    r"insertion sequence|IS[0-9]|excisionase|prophage|tape measure|holin|"
    r"lysin|antirepressor|baseplate", re.I)

# STRING channel column indices in *.protein.links.detailed.*.txt (SPACE-delimited)
# 0 protein1  1 protein2  2 neighborhood  3 fusion  4 cooccurence
# 5 coexpression  6 experimental  7 database  8 textmining  9 combined_score
CH = {"neighborhood": 2, "fusion": 3, "cooccurence": 4, "coexpression": 5,
      "experimental": 6, "database": 7, "textmining": 8, "combined_score": 9}
POSITIONAL = ["neighborhood", "fusion", "cooccurence"]
PRIMARY = ["experimental", "database", "textmining"]
SENSITIVITY = ["experimental", "database", "textmining", "coexpression"]

# Alias sources that carry a gene SYMBOL (preferred join key: symbols are
# conserved across strains, whereas locus tags differ between our assemblies
# and STRING's reference strain).
SYMBOL_SOURCES = {"RefSeq_gene", "UniProt_GN_Name", "KEGG_NAME",
                  "BLAST_UniProt_GN_Name", "Ensembl_gene", "Ensembl_EntrezGene"}
LOCUS_SOURCES = {"RefSeq_locus", "UniProt_GN_OrderedLocusNames",
                 "RefSeq_synonym", "Ensembl_locus"}

STRING_PRIOR = 0.041  # STRING's random-expectation prior, per their docs


def combine_channels(scores):
    """STRING's own probabilistic channel combination, applied to a SUBSET of
    channels. Each raw score is 0-1000. Removes the random prior from each
    channel, combines as independent evidence (probabilistic OR), re-adds the
    prior. Returns 0-1000."""
    prod = 1.0
    for raw in scores:
        s = raw / 1000.0
        s_adj = (s - STRING_PRIOR) / (1 - STRING_PRIOR)
        if s_adj < 0:
            s_adj = 0.0
        prod *= (1 - s_adj)
    combined = 1 - prod
    combined = combined + STRING_PRIOR * (1 - combined)
    return round(combined * 1000, 1)


def load_aliases(taxid):
    """Return (symbol_map, locus_map): lowercased alias -> set of STRING ids."""
    path = os.path.join(SDIR, f"{taxid}.protein.aliases.{STRING_VERSION}.txt.gz")
    symbol_map, locus_map = {}, {}
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            sid, alias, source = parts[0], parts[1], parts[2]
            if not alias:
                continue
            key = alias.lower()
            if source in SYMBOL_SOURCES:
                symbol_map.setdefault(key, set()).add(sid)
            elif source in LOCUS_SOURCES:
                locus_map.setdefault(key, set()).add(sid)
    return symbol_map, locus_map


def load_links(taxid, wanted_ids):
    """Stream the links file, keeping only interactions where BOTH partners are
    in wanted_ids. Returns dict frozenset({id1,id2}) -> list of raw channel scores."""
    path = os.path.join(SDIR, f"{taxid}.protein.links.detailed.{STRING_VERSION}.txt.gz")
    links = {}
    with gzip.open(path, "rt") as f:
        header = f.readline()  # space-delimited header
        for line in f:
            parts = line.split()
            if len(parts) < 10:
                continue
            p1, p2 = parts[0], parts[1]
            if p1 not in wanted_ids or p2 not in wanted_ids:
                continue
            links[frozenset((p1, p2))] = [int(x) for x in parts[2:10]]
    return links


def resolve(gene, symbol_map, locus_map, ortho):
    """Map one GFF gene record to STRING ids.

    Order of preference:
      1. BLAST-based ortholog map (build_ortholog_map.py) -- sequence homology,
         works even when STRING's reference strain shares no ID namespace with
         our assembly (the V. cholerae / S. pyogenes / B. fragilis case).
      2. Gene-symbol alias join -- fine where STRING carries symbols.
      3. Locus-tag alias join -- last resort; only works if STRING's reference
         strain happens to use the same locus tags as our assembly.
    """
    sym = gene.get("symbol")
    if sym and ortho:
        sid = ortho.get(sym)
        if sid:
            return {sid}, "ortholog"
    if sym:
        ids = symbol_map.get(sym.lower())
        if ids:
            return ids, "symbol"
    for key in (gene.get("locus_tag"), sym):
        if key:
            ids = locus_map.get(key.lower())
            if ids:
                return ids, "locus_tag"
    return None, None


def main():
    with open(os.path.join(DATA_DIR, "divergent_disruption_results.json")) as f:
        scan = json.load(f)

    taxids = {}
    with open(os.path.join(SDIR, "taxids.tsv")) as f:
        for line in f:
            name, taxid = line.strip().split("\t")
            taxids[name] = taxid

    # BLAST-based ortholog map, if build_ortholog_map.py has been run
    ortho_path = os.path.join(DATA_DIR, "ortholog_map.json")
    ortho_all = {}
    if os.path.exists(ortho_path):
        with open(ortho_path) as f:
            ortho_all = json.load(f).get("map", {})
        print(f"using BLAST ortholog map for {len(ortho_all)} species")
    else:
        print("WARNING: ortholog_map.json not found -- symbol/locus join only")

    rows = []
    per_species_stats = {}

    for species in sorted(taxids):
        sp_path = os.path.join(GDIR, species)
        if not os.path.isdir(sp_path):
            continue
        taxid = taxids[species]

        # --- reference genome & its divergent pairs (same definition as Parts 1/3)
        ref_acc = fdd.REFERENCE_STRAIN.get(species)
        accs = sorted(os.listdir(sp_path))
        if ref_acc not in accs:
            ref_acc = accs[0]
        ref_genes, _ = fdd.load_gff(os.path.join(sp_path, ref_acc, "genomic.gff"))
        ref_lookup = fdd.gene_lookup(ref_genes)
        pairs = fdd.divergent_pairs(ref_genes, max_len=500)

        # --- disruption status per pair, from the Part 1/3 scan
        disrupt = {}
        for h in scan.get(species, {}).get("hits", []):
            if not (MIN_SIZE <= h["expanded_length"] <= MAX_SIZE):
                continue
            key = tuple(sorted(h["pair"]))
            d = disrupt.setdefault(key, {"n": 0, "max_len": 0, "mobile": False})
            d["n"] += 1
            d["max_len"] = max(d["max_len"], h["expanded_length"])
            if re.search(MOBILE_KW, " ".join(p or "" for _, p in h["intervening_products"])):
                d["mobile"] = True

        # --- STRING id resolution
        symbol_map, locus_map = load_aliases(taxid)
        ortho = ortho_all.get(species, {})
        resolved, how = {}, {}
        for (a, b) in pairs:
            for sym in (a, b):
                if sym in resolved:
                    continue
                g = ref_lookup.get(sym, [None])[0]
                if g is None:
                    resolved[sym] = None
                    continue
                ids, method = resolve(g, symbol_map, locus_map, ortho)
                resolved[sym] = ids
                how[sym] = method

        wanted = set()
        for ids in resolved.values():
            if ids:
                wanted |= ids
        links = load_links(taxid, wanted)

        n_mapped_pairs = 0
        for (a, b), info in pairs.items():
            ids_a, ids_b = resolved.get(a), resolved.get(b)
            mapped = bool(ids_a and ids_b)
            row = {
                "species": species, "taxid": taxid,
                "geneA": a, "geneB": b, "spacer": info["length"],
                "mapped": mapped,
                "map_method_A": how.get(a), "map_method_B": how.get(b),
                "disrupted": False,
                "n_strains_disrupted": 0, "max_expanded_len": None, "mobile_kw": False,
            }
            key = tuple(sorted((a, b)))
            if key in disrupt:
                row["disrupted"] = True
                row["n_strains_disrupted"] = disrupt[key]["n"]
                row["max_expanded_len"] = disrupt[key]["max_len"]
                row["mobile_kw"] = disrupt[key]["mobile"]

            if not mapped:
                # unmappable: leave scores as None. Critically NOT zero -- "no
                # annotation" and "annotated as uncoupled" are different states
                # and conflating them would bias the selection test.
                for ch in CH:
                    row[ch] = None
                row["score_primary"] = None
                row["score_with_coexpr"] = None
                rows.append(row)
                continue

            n_mapped_pairs += 1
            # best-scoring combination if a symbol is ambiguous across ids
            best = None
            for ia in ids_a:
                for ib in ids_b:
                    if ia == ib:
                        continue
                    ch_scores = links.get(frozenset((ia, ib)))
                    if ch_scores is None:
                        continue
                    prim = combine_channels([ch_scores[CH[c] - 2] for c in PRIMARY])
                    if best is None or prim > best[0]:
                        best = (prim, ch_scores)
            if best is None:
                # both genes mapped, but STRING lists no interaction between
                # them at all -> genuine zero evidence in every channel
                for ch in CH:
                    row[ch] = 0
                row["score_primary"] = combine_channels([0, 0, 0])
                row["score_with_coexpr"] = combine_channels([0, 0, 0, 0])
            else:
                prim, ch_scores = best
                for ch, idx in CH.items():
                    row[ch] = ch_scores[idx - 2]
                row["score_primary"] = prim
                row["score_with_coexpr"] = combine_channels(
                    [ch_scores[CH[c] - 2] for c in SENSITIVITY])
            rows.append(row)

        method_counts = {}
        for m in how.values():
            method_counts[m or "unmapped"] = method_counts.get(m or "unmapped", 0) + 1
        per_species_stats[species] = {
            "taxid": taxid, "n_pairs": len(pairs), "n_mapped": n_mapped_pairs,
            "map_rate": round(n_mapped_pairs / len(pairs), 3) if pairs else None,
            "gene_resolution_methods": method_counts,
        }
        print(f"{species:<14} taxid={taxid:<7} pairs={len(pairs):<5} "
              f"STRING-mapped={n_mapped_pairs:<5} ({per_species_stats[species]['map_rate']:.0%})  "
              f"methods={method_counts}")

    out = {"per_species": per_species_stats, "pairs": rows,
           "channels": {"excluded_positional": POSITIONAL,
                        "primary": PRIMARY, "sensitivity": SENSITIVITY}}
    with open(os.path.join(DATA_DIR, "coupling_scores.json"), "w") as f:
        json.dump(out, f, indent=2)

    # analysis-ready flat table
    cols = ["species", "geneA", "geneB", "spacer", "mapped", "score_primary",
            "score_with_coexpr", "neighborhood", "fusion", "cooccurence",
            "coexpression", "experimental", "database", "textmining",
            "combined_score", "disrupted", "n_strains_disrupted",
            "max_expanded_len", "mobile_kw"]
    with open(os.path.join(DATA_DIR, "coupling_scores.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")

    total = len(rows)
    mapped = sum(1 for r in rows if r["mapped"])
    print(f"\nTOTAL: {total} divergent pairs, {mapped} STRING-mapped ({mapped/total:.0%})")
    print(f"wrote coupling_scores.json / coupling_scores.tsv")


if __name__ == "__main__":
    main()
