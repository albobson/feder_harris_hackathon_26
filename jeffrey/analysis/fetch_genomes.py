import json, os, time, urllib.request

# Repo-relative paths: species list lives in ../data, downloaded genomes go there too.
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(os.path.join(HERE, "..", "data"))
GDIR = os.path.join(DATA_DIR, "genomes")
N_PER_SPECIES = 4  # genomes (strains) to pull per species

species = []
with open(os.path.join(DATA_DIR, "species_list.txt")) as f:
    for line in f:
        taxid, name = line.strip().split("\t")
        species.append((taxid, name))

def fetch_json(url):
    req = urllib.request.Request(url, headers={"User-Agent": "hackathon-script"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read())

manifest = []
for taxid, name in species:
    # NCBI Datasets v2 REST API: report of RefSeq/GenBank assemblies for this taxon,
    # restricted to finished (Complete Genome) assemblies so contigs/scaffolds aren't
    # a confound when we later compare intergenic distances across strains.
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{taxid}/dataset_report?filters.assembly_level=complete_genome&page_size=40"
    try:
        data = fetch_json(url)
    except Exception as e:
        print(f"ERROR querying {name}: {e}")
        continue
    reports = data.get("reports", [])
    picked = []
    seen_strains = set()
    for rep in reports:
        acc = rep.get("accession")
        # Prefer the RefSeq (GCF_) accession when a GenBank/RefSeq pair exists,
        # since RefSeq annotation (gene symbols) is more consistent across strains.
        paired = rep.get("assembly_info", {}).get("paired_assembly", {}).get("accession")
        use_acc = paired if paired and paired.startswith("GCF") else acc
        strain = rep.get("organism", {}).get("infraspecific_names", {}).get("strain", use_acc)
        if strain in seen_strains:
            continue  # skip duplicate re-submissions of the same strain
        seen_strains.add(strain)
        picked.append((use_acc, strain))
        if len(picked) >= N_PER_SPECIES:
            break
    print(f"{name} ({taxid}): found {len(reports)} complete genomes, picked {len(picked)}")
    for acc, strain in picked:
        manifest.append((name, taxid, acc, strain))
    time.sleep(0.4)  # be polite to the NCBI API (no API key used)

with open(os.path.join(DATA_DIR, "genome_manifest.tsv"), "w") as f:
    for row in manifest:
        f.write("\t".join(row) + "\n")

print(f"\nTotal genomes to download: {len(manifest)}")
