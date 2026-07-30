"""Part 2 of the survey: instead of starting from "any short-spacer divergent
pair" and asking whether it's disrupted, start from divergent pairs that are
PLAUSIBLY FUNCTIONALLY COUPLED (a transcriptional regulator divergently
transcribed from the operon it likely controls -- the same general
architecture as torS/torT, araC/araBAD, acrR/acrAB, soxR/soxS, etc., where the
regulator autoregulates itself and/or its target via a shared promoter/operator)
and ask whether those specific pairs are ever disrupted by an insertion.

Reuses the GFF-parsing/divergent-pair/circular-span helpers from
find_divergent_disruptions.py (imported, not duplicated) plus the disruption
results it already computed and saved to divergent_disruption_results.json.
"""
import os, re, json
import find_divergent_disruptions as fdd

DATA_DIR = fdd.DATA_DIR
GDIR = fdd.GDIR

# General heuristic for "this divergent pair is plausibly functionally coupled":
# one side's product annotation names a regulator family/role. This is the same
# broad architectural pattern that makes torS/torT, araC/araBAD, acrR/acrAB,
# soxR/soxS, and mepR/mepA functionally coupled -- a regulator sharing its
# promoter/operator with the thing it regulates -- generalized so it doesn't
# require hand-curating a gene-symbol list for every species.
REGULATOR_KW = re.compile(
    r"transcriptional (regulator|repressor|activator)|regulatory protein|"
    r"(TetR|LysR|MarR|AraC|GntR|IclR|LacI|ArsR|Crp|Fnr|MerR|DeoR|Fur|PadR|XRE)[- ]family|"
    r"response regulator|sigma factor|anti-sigma|two-component system|"
    r"repressor protein|activator protein",
    re.I,
)

# Curated positive examples with precise, literature-verified gene symbols and
# functional relationships (filled in from a verification pass -- see SUMMARY.md
# Part 4 for sourcing/caveats on each). Checked directly against every species'
# genomes below, not just the species where the relationship was first described,
# since some of these pairs are conserved across many Gammaproteobacteria.
CURATED_PAIRS = [
    ("torS", "torT", "same TMAO-sensing two-component system (Carey et al. 2019)"),
    ("metE", "metR", "MetR (LysR-family) activates metE transcription (methionine biosynthesis)"),
    ("araC", "araB", "AraC autoregulates and activates araBAD (arabinose catabolism) from a shared control region"),
    ("soxR", "soxS", "SoxR (redox sensor) directly activates soxS (oxidative-stress regulon amplifier)"),
    ("acrR", "acrA", "AcrR represses acrAB(-tolC) multidrug efflux operon from a shared promoter"),
    ("mepR", "mepA", "MepR (MarR-family) represses mepA multidrug efflux pump (S. aureus, chromosomal)"),
    ("cpxP", "cpxR", "CpxP (periplasmic accessory inhibitor) divergently transcribed from cpxRA operon it helps regulate, conserved into V. cholerae"),
]

def load_species_reference(species_dir):
    sp_path = os.path.join(GDIR, species_dir)
    ref_acc = fdd.REFERENCE_STRAIN.get(species_dir)
    accs = sorted(os.listdir(sp_path))
    if ref_acc not in accs:
        ref_acc = accs[0]
    gff = os.path.join(sp_path, ref_acc, "genomic.gff")
    genes, seq_lengths = fdd.load_gff(gff)
    return ref_acc, genes, seq_lengths, accs

def check_pair_across_strains(species_dir, symA, symB, sp_path, accs, ref_acc, baseline_len):
    """Same disruption logic as the main scan, applied to one specific gene-symbol
    pair across every non-reference strain of a species."""
    hits = []
    for acc in accs:
        if acc == ref_acc:
            continue
        gff = os.path.join(sp_path, acc, "genomic.gff")
        if not os.path.exists(gff):
            continue
        genes, seq_lengths = fdd.load_gff(gff)
        lookup = fdd.gene_lookup(genes)
        span_info = fdd.find_span(lookup, symA, symB, seq_lengths)
        if span_info is None:
            continue
        if span_info["span"] > 5 * max(baseline_len, 50) and span_info["span"] > 800 and span_info["span"] <= 200_000:
            prods = fdd.intervening_products(genes, span_info["seqid"],
                                              span_info["left"]["end"], span_info["right"]["start"],
                                              span_info["length"], span_info["wrapped"])
            hits.append({"genome": acc, "ref_length": baseline_len, "expanded_length": span_info["span"],
                         "n_intervening_genes": len(prods), "intervening_products": prods})
    return hits

def main():
    # Load the already-computed disruption results from the main scan so the
    # heuristic pass can just ask "was this pair already found disrupted?"
    # instead of recomputing the whole cross-strain scan.
    with open(os.path.join(DATA_DIR, "divergent_disruption_results.json")) as f:
        existing_results = json.load(f)

    heuristic_report = {}
    curated_report = {}

    for species_dir in sorted(os.listdir(GDIR)):
        sp_path = os.path.join(GDIR, species_dir)
        if not os.path.isdir(sp_path):
            continue
        ref_acc, ref_genes, ref_seqlens, accs = load_species_reference(species_dir)
        ref_lookup = fdd.gene_lookup(ref_genes)
        ref_pairs = fdd.divergent_pairs(ref_genes, max_len=500)

        # --- heuristic regulator-headed pairs ---
        regulator_pairs = []
        for (symA, symB), info in ref_pairs.items():
            prodA = ref_lookup[symA][0]["product"] if symA in ref_lookup else None
            prodB = ref_lookup[symB][0]["product"] if symB in ref_lookup else None
            is_reg_A = bool(prodA and REGULATOR_KW.search(prodA))
            is_reg_B = bool(prodB and REGULATOR_KW.search(prodB))
            if is_reg_A or is_reg_B:
                regulator_pairs.append({
                    "pair": (symA, symB), "spacer": info["length"],
                    "regulator_gene": symA if is_reg_A else symB,
                    "regulator_product": prodA if is_reg_A else prodB,
                    "partner_product": prodB if is_reg_A else prodA,
                })

        existing_hits_by_pair = {}
        for h in existing_results.get(species_dir, {}).get("hits", []):
            key = tuple(sorted(h["pair"]))
            existing_hits_by_pair.setdefault(key, []).append(h)

        annotated = []
        for rp in regulator_pairs:
            key = tuple(sorted(rp["pair"]))
            rp["disrupted_in"] = existing_hits_by_pair.get(key, [])
            annotated.append(rp)
        heuristic_report[species_dir] = {
            "reference": ref_acc,
            "n_regulator_headed_divergent_pairs": len(annotated),
            "pairs": annotated,
        }

        # --- curated pairs ---
        curated_hits = []
        for symA, symB, relationship in CURATED_PAIRS:
            span_info = fdd.find_span(ref_lookup, symA, symB, ref_seqlens)
            if span_info is None:
                continue  # this pair isn't present / not divergently oriented in this species
            baseline_len = span_info["span"]
            strain_hits = check_pair_across_strains(species_dir, symA, symB, sp_path, accs, ref_acc, baseline_len)
            curated_hits.append({
                "pair": (symA, symB), "relationship": relationship,
                "baseline_spacer": baseline_len, "disrupted_in": strain_hits,
            })
        if curated_hits:
            curated_report[species_dir] = {"reference": ref_acc, "pairs": curated_hits}

    with open(os.path.join(DATA_DIR, "functional_pair_results.json"), "w") as f:
        json.dump({"heuristic": heuristic_report, "curated": curated_report}, f, indent=2)

    print("=== Heuristic regulator-headed divergent pairs (per species) ===")
    for sp, d in heuristic_report.items():
        n_disrupted = sum(1 for p in d["pairs"] if p["disrupted_in"])
        print(f"{sp}: {d['n_regulator_headed_divergent_pairs']} regulator-headed divergent pairs, "
              f"{n_disrupted} disrupted in >=1 sampled strain")

    print("\n=== Curated functionally-coupled pairs (per species present) ===")
    for sp, d in curated_report.items():
        for p in d["pairs"]:
            status = f"DISRUPTED in {len(p['disrupted_in'])} strain(s)" if p["disrupted_in"] else "not disrupted in this sample"
            print(f"{sp}: {p['pair']} (baseline {p['baseline_spacer']}bp) -- {status}")

if __name__ == "__main__":
    main()
