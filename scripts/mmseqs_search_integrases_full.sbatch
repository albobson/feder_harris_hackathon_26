#!/bin/bash
#SBATCH --job-name=mmseqs_integrases_full
#SBATCH --partition=campus-new
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=08:00:00
#SBATCH --output=logs/mmseqs_integrases_full_%j.out
#SBATCH --error=logs/mmseqs_integrases_full_%j.err

# Full-DB CPU search: 10 phage integrases -> OMG clustered90 reps.
# See mmseqs_search_torT_torS_full.sbatch for split/mem rationale.
# Integrase --max-seqs bumped from 20000 (pilot) to 500000: several
# integrase queries came close to the pilot cap (e.g. HK620 6277,
# Sf6 6325), and tyrosine recombinases are extremely broad in
# metagenomes.

set -euo pipefail

REPO=/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/feder_harris_hackathon_26
DB=/fh/working/srivatsan_s/databases/mmseqs2/omg_cds/omg_cds_clustered90_seqs

cd "$REPO"
mkdir -p results tmp logs

ml MMseqs2/13-45111-gompi-2021b
mmseqs version

mmseqs createdb queries/phage_integrases.faa queries/phage_integrases_db

mmseqs search \
    queries/phage_integrases_db \
    "$DB" \
    results/integrases_vs_omgfull \
    tmp/integrases_full \
    -s 7.5 \
    -e 1e-10 \
    --max-seqs 500000 \
    --split-memory-limit 200G \
    --threads "$SLURM_CPUS_PER_TASK"

mmseqs convertalis \
    queries/phage_integrases_db \
    "$DB" \
    results/integrases_vs_omgfull \
    results/integrases_vs_omgfull.m8 \
    --format-output "query,target,theader,pident,alnlen,evalue,bits,qcov,tcov"

echo "done. hits: $(wc -l < results/integrases_vs_omgfull.m8)"
