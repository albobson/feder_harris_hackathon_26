import json, os, random, time, urllib.request

# Repo-relative paths: species list lives in ../data, downloaded genomes go there too.
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.normpath(os.path.join(HERE, "..", "data"))
GDIR = os.path.join(DATA_DIR, "genomes")
N_PER_SPECIES = 20     # genomes (strains) to pull per species
POOL_PAGE_SIZE = 1000  # size of the candidate pool to sample from, per species
RANDOM_SEED = 20260730 # fixed seed so the "random" sample is reproducible across reruns

species = []
with open(os.path.join(DATA_DIR, "species_list.txt")) as f:
    for line in f:
        taxid, name = line.strip().split("\t")
        species.append((taxid, name))

def fetch_json(url):
    req = urllib.request.Request(url, headers={"User-Agent": "hackathon-script"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read())

rng = random.Random(RANDOM_SEED)
manifest = []
for taxid, name in species:
    # NCBI Datasets v2 REST API: report of RefSeq/GenBank assemblies for this taxon,
    # restricted to finished (Complete Genome) assemblies so contigs/scaffolds aren't
    # a confound when we later compare intergenic distances across strains. Pull a
    # large pool (up to POOL_PAGE_SIZE) so we can randomly sample strains from it,
    # rather than just taking whatever order the API returns them in -- the API's
    # default order isn't documented as random, and could cluster genomes submitted
    # together by the same study/sequencing center.
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{taxid}/dataset_report?filters.assembly_level=complete_genome&page_size={POOL_PAGE_SIZE}"
    try:
        data = fetch_json(url)
    except Exception as e:
        print(f"ERROR querying {name}: {e}")
        continue
    reports = data.get("reports", [])

    # Dedup by strain (a strain can have both a GenBank and paired RefSeq entry;
    # keep one accession per distinct strain), preferring the RefSeq (GCF_)
    # accession since its annotation (gene symbols) is more consistent across strains.
    by_strain = {}
    for rep in reports:
        acc = rep.get("accession")
        paired = rep.get("assembly_info", {}).get("paired_assembly", {}).get("accession")
        use_acc = paired if paired and paired.startswith("GCF") else acc
        strain = rep.get("organism", {}).get("infraspecific_names", {}).get("strain", use_acc)
        by_strain.setdefault(strain, use_acc)

    pool = list(by_strain.items())  # [(strain, accession), ...]
    rng.shuffle(pool)
    picked = pool[:N_PER_SPECIES]
    print(f"{name} ({taxid}): found {len(reports)} complete genome records, "
          f"{len(pool)} distinct strains, randomly picked {len(picked)}")
    for strain, acc in picked:
        manifest.append((name, taxid, acc, strain))
    time.sleep(0.4)  # be polite to the NCBI API (no API key used)

with open(os.path.join(DATA_DIR, "genome_manifest.tsv"), "w") as f:
    for row in manifest:
        f.write("\t".join(row) + "\n")

print(f"\nTotal genomes to download: {len(manifest)}")
