"""Re-run the Part 5 selection test using the validated STRING coupling score
instead of the regulator-keyword heuristic that made the sign of the effect flip.

Only pairs that MAP to STRING are used. Unmappable pairs are dropped rather than
scored zero -- treating "no annotation" as "not coupled" is what would bias the
result.
"""
import csv, math, os
from math import comb
import find_divergent_disruptions as fdd

DATA_DIR = fdd.DATA_DIR
COUPLED_THRESHOLD = 400.0   # STRING's conventional "medium confidence" cutoff


def fisher_two_sided(a, b, c, d):
    r1, r2, c1 = a + b, c + d, a + c
    N = a + b + c + d
    denom = comb(N, c1)
    def P(k):
        return comb(r1, k) * comb(r2, c1 - k) / denom
    p_obs = P(a)
    lo, hi = max(0, c1 - r2), min(r1, c1)
    return sum(P(k) for k in range(lo, hi + 1) if P(k) <= p_obs * (1 + 1e-9))


def or_ci(a, b, c, d):
    # Haldane-Anscombe correction when any cell is zero
    if 0 in (a, b, c, d):
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    odds = (a * d) / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return odds, math.exp(math.log(odds) - 1.96 * se), math.exp(math.log(odds) + 1.96 * se)


def median(xs):
    if not xs:
        return None
    s = sorted(xs); m = len(s) // 2
    return s[m] if len(s) % 2 else (s[m - 1] + s[m]) / 2


rows = []
with open(os.path.join(DATA_DIR, "coupling_scores.tsv")) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["score_primary"] in ("", None):
            continue                      # unmappable -> excluded, NOT scored zero
        rows.append({
            "species": r["species"],
            "score": float(r["score_primary"]),
            "spacer": int(r["spacer"]),
            "disrupted": r["disrupted"] == "True",
        })

coupled = [r for r in rows if r["score"] >= COUPLED_THRESHOLD]
uncoupled = [r for r in rows if r["score"] < COUPLED_THRESHOLD]
a = sum(1 for r in coupled if r["disrupted"])       # coupled & disrupted
b = len(coupled) - a
c = sum(1 for r in uncoupled if r["disrupted"])     # uncoupled & disrupted
d = len(uncoupled) - c

print("=" * 78)
print("SELECTION TEST with validated STRING coupling score (AUC 0.872 vs hand labels)")
print("=" * 78)
print(f"mapped pairs analysed        : {len(rows)} of 4099 ({len(rows)/4099:.0%})")
print(f"coupling threshold           : score_primary >= {COUPLED_THRESHOLD:.0f} (STRING medium confidence)")
print()
print(f"{'':<26}{'disrupted':>10}{'not':>10}{'rate':>9}")
print(f"{'functionally COUPLED':<26}{a:>10}{b:>10}{a/(a+b)*100:>8.2f}%")
print(f"{'not coupled':<26}{c:>10}{d:>10}{c/(c+d)*100:>8.2f}%")
odds, lo, hi = or_ci(a, b, c, d)
p = fisher_two_sided(a, b, c, d)
print()
print(f"odds ratio = {odds:.2f}   95% CI [{lo:.2f}, {hi:.2f}]   Fisher exact p = {p:.3f}")
direction = "ENRICHED" if odds > 1 else "DEPLETED"
sig = "significant" if p < 0.05 else "NOT significant"
print(f"-> disruption is {direction} at functionally-coupled pairs; {sig} at alpha=0.05")

print()
print("-" * 78)
print("CONFOUND CHECK: is spacer length different between the two groups?")
print("-" * 78)
print(f"median spacer, coupled     = {median([r['spacer'] for r in coupled])} bp")
print(f"median spacer, not coupled = {median([r['spacer'] for r in uncoupled])} bp")
print("(if these differ substantially, spacer length -- which plausibly drives")
print(" insertion independent of function -- needs to be controlled for)")

print()
print("-" * 78)
print("PER-SPECIES BREAKDOWN (species differ ~100x in disruption rate)")
print("-" * 78)
print(f"{'species':<15}{'coupled d/n':>14}{'uncoupled d/n':>16}")
for sp in sorted(set(r["species"] for r in rows)):
    ca = sum(1 for r in coupled if r["species"] == sp and r["disrupted"])
    cn = sum(1 for r in coupled if r["species"] == sp)
    ua = sum(1 for r in uncoupled if r["species"] == sp and r["disrupted"])
    un = sum(1 for r in uncoupled if r["species"] == sp)
    print(f"{sp:<15}{f'{ca}/{cn}':>14}{f'{ua}/{un}':>16}")

print()
print("NOTE: with only ~26 disruption events total this test remains severely")
print("underpowered -- see Part 5. The annotation fix makes the ESTIMATE")
print("trustworthy; it does not make the SAMPLE adequate.")
