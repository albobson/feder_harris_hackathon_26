"""Design/power analysis: how many genomes would we need to say something
provisional about SELECTION for/against insertion elements landing between
divergent, functionally-related gene pairs?

Uses the empirical base rates from our own scans (Parts 3/4) rather than
guessed numbers. Everything here is an order-of-magnitude estimate: the
base rate itself is estimated from only ~26 disruption events, so its own
confidence interval is wide and propagates into every N below.
"""
import json, math, os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(os.path.join(HERE, "..", "data"))

MIN_SIZE, MAX_SIZE = 800, 200_000  # discrete-insertion size window used throughout

with open(os.path.join(DATA_DIR, "divergent_disruption_results.json")) as f:
    scan = json.load(f)
with open(os.path.join(DATA_DIR, "functional_pair_results.json")) as f:
    fp = json.load(f)["heuristic"]

# ---------------------------------------------------------------- base rates
n_pairs_total = 0        # divergent pairs monitored (reference-strain pairs)
n_pairs_reg = 0          # of those, "regulator-headed" (functional-coupling proxy)
disrupted_pairs = set()  # unique (species, pair) disrupted in >=1 strain
disrupted_strain_level = 0  # strain-level hits, i.e. counting each strain separately
n_strains = 0

for sp, d in scan.items():
    n_pairs_total += d["n_ref_pairs"]
    hits = [h for h in d["hits"] if MIN_SIZE <= h["expanded_length"] <= MAX_SIZE]
    disrupted_strain_level += len(hits)
    for h in hits:
        disrupted_pairs.add((sp, tuple(sorted(h["pair"]))))
    # count sampled strains actually present on disk for this species
    sp_dir = os.path.join(DATA_DIR, "genomes", sp)
    if os.path.isdir(sp_dir):
        n_strains += len([a for a in os.listdir(sp_dir)
                          if os.path.exists(os.path.join(sp_dir, a, "genomic.gff"))])

reg_disrupted_pairs = set()
for sp, d in fp.items():
    n_pairs_reg += d["n_regulator_headed_divergent_pairs"]
    for p in d["pairs"]:
        good = [h for h in p["disrupted_in"] if MIN_SIZE <= h["expanded_length"] <= MAX_SIZE]
        if good:
            reg_disrupted_pairs.add((sp, tuple(sorted(p["pair"]))))

n_species = len(scan)
strains_per_species = n_strains / n_species
# non-reference strains are the ones that can show a disruption vs the baseline
comparisons_per_pair = strains_per_species - 1

# Literature-verification survival rate: of the regulator-headed pairs that
# showed disruption and were checked by hand, how many were genuinely
# functionally coupled? rocD/rocR + yyaT/yybA held up; btsS/mlrA, purT/ybfI,
# fabR/sthA did not (3 more were coordinate artifacts, already size-filtered).
VERIFIED_TRUE, VERIFIED_CHECKED = 2, 5
verify_survival = VERIFIED_TRUE / VERIFIED_CHECKED

p_any = len(disrupted_pairs) / n_pairs_total
p_reg = len(reg_disrupted_pairs) / n_pairs_reg
p_reg_verified = p_reg * verify_survival
n_pairs_nonreg = n_pairs_total - n_pairs_reg
p_nonreg = (len(disrupted_pairs) - len(reg_disrupted_pairs)) / n_pairs_nonreg

print("=" * 72)
print("EMPIRICAL BASE RATES (from our 180-genome scan)")
print("=" * 72)
print(f"species sampled                         : {n_species}")
print(f"genomes with usable annotation          : {n_strains} ({strains_per_species:.1f}/species)")
print(f"divergent pairs monitored (all)         : {n_pairs_total}")
print(f"  ...regulator-headed (coupling proxy)  : {n_pairs_reg}")
print(f"  ...non-regulator-headed               : {n_pairs_nonreg}")
print()
print(f"strain-level disruption hits            : {disrupted_strain_level}")
print(f"UNIQUE pairs disrupted in >=1 strain    : {len(disrupted_pairs)}")
print(f"  -> pseudo-replication collapse factor : {disrupted_strain_level/len(disrupted_pairs):.1f}x")
print()
print(f"P(disrupted | any divergent pair)       : {p_any:.5f}  ({p_any*100:.3f}%)")
print(f"P(disrupted | regulator-headed)         : {p_reg:.5f}  ({p_reg*100:.3f}%)")
print(f"P(disrupted | non-regulator-headed)     : {p_nonreg:.5f}  ({p_nonreg*100:.3f}%)")
print(f"P(disrupted | VERIFIED coupled)         : {p_reg_verified:.5f}  ({p_reg_verified*100:.3f}%)")
print(f"  (applying {VERIFIED_TRUE}/{VERIFIED_CHECKED} literature-verification survival rate)")

# --------------------------------------------------- two-proportion power calc
Z_ALPHA = 1.959964  # two-sided alpha = 0.05
Z_BETA = 0.841621   # power = 0.80

def n_per_group(p0, p1):
    """Pairs needed per group to distinguish proportions p0 vs p1."""
    if p0 == p1:
        return float("inf")
    pbar = (p0 + p1) / 2
    num = (Z_ALPHA * math.sqrt(2 * pbar * (1 - pbar))
           + Z_BETA * math.sqrt(p0 * (1 - p0) + p1 * (1 - p1))) ** 2
    return num / (p0 - p1) ** 2

print()
print("=" * 72)
print("TIER 3: COMPARATIVE SELECTION TEST (coupled vs non-coupled pairs)")
print("=" * 72)
print("How many divergent pairs must be monitored PER GROUP to detect a given")
print("effect, at 80% power / alpha=0.05, holding strains-per-species at ~20?")
print()
print(f"{'effect on coupled pairs':<34} {'pairs/group':>12} {'species':>9} {'genomes':>9}")
print("-" * 72)
pairs_per_species_reg = n_pairs_reg / n_species
for label, p1 in [
    ("5x depletion (strong purifying)", p_nonreg / 5),
    ("2x depletion (moderate)",         p_nonreg / 2),
    ("1.5x depletion (weak)",           p_nonreg / 1.5),
    ("2x enrichment (adaptive/hotspot)", p_nonreg * 2),
    ("5x enrichment (strong)",          p_nonreg * 5),
]:
    n = n_per_group(p_nonreg, p1)
    species_needed = n / pairs_per_species_reg
    genomes = species_needed * strains_per_species
    print(f"{label:<34} {n:>12,.0f} {species_needed:>9,.0f} {genomes:>9,.0f}")

print()
print("Same, but if we instead deepen sampling to 100 strains/species (which")
print("raises per-pair detection probability ~5x, to the extent insertions are")
print("polymorphic within a species rather than fixed):")
print()
scale = 5.0
p_nonreg_deep = min(p_nonreg * scale, 0.5)
print(f"{'effect on coupled pairs':<34} {'pairs/group':>12} {'species':>9} {'genomes':>9}")
print("-" * 72)
for label, div in [("5x depletion (strong purifying)", 5),
                   ("2x depletion (moderate)", 2),
                   ("2x enrichment (adaptive/hotspot)", 0.5)]:
    p1 = p_nonreg_deep / div
    n = n_per_group(p_nonreg_deep, p1)
    species_needed = n / pairs_per_species_reg
    genomes = species_needed * 100
    print(f"{label:<34} {n:>12,.0f} {species_needed:>9,.0f} {genomes:>9,.0f}")

# ------------------------------------------- Tier 2: just estimating the rate
print()
print("=" * 72)
print("TIER 2: PRECISE RATE ESTIMATE (no comparison, just 'how often?')")
print("=" * 72)
print("Poisson relative precision needs k events for +/-1/sqrt(k) rel. SE:")
print()
# current yield of VERIFIED coupled-pair disruption events
current_verified_events = len(reg_disrupted_pairs) * verify_survival
print(f"verified coupled-pair events in hand    : ~{current_verified_events:.1f}")
print()
print(f"{'target precision':<24} {'events needed':>14} {'scale-up':>10} {'genomes':>10}")
print("-" * 72)
for label, k in [("+/-50% (very rough)", 4), ("+/-30% (rough)", 11),
                 ("+/-20% (usable)", 25), ("+/-10% (tight)", 100)]:
    scale_up = k / current_verified_events
    genomes = scale_up * n_strains
    print(f"{label:<24} {k:>14} {scale_up:>9.0f}x {genomes:>10,.0f}")

print()
print("=" * 72)
print("KEY CAVEAT: EFFECTIVE SAMPLE SIZE != GENOME COUNT")
print("=" * 72)
print(f"Our {disrupted_strain_level} strain-level hits collapse to {len(disrupted_pairs)} unique pairs")
print(f"({disrupted_strain_level/len(disrupted_pairs):.1f}x pseudo-replication). Strains sharing an insertion by")
print("common descent are ONE evolutionary event, not many. Every N above is")
print("therefore a floor: it assumes independent events. Correcting for")
print("phylogeny (counting insertions per independent lineage on a core-genome")
print("tree) will inflate required genome counts further, plausibly 2-5x.")
