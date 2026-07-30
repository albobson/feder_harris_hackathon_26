#!/bin/bash
#SBATCH --job-name=mmseqs_torT_torS_full
#SBATCH --partition=campus-new
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=08:00:00
#SBATCH --output=logs/mmseqs_torT_torS_full_%j.out
#SBATCH --error=logs/mmseqs_torT_torS_full_%j.err

# Full-DB CPU search: MG1655 torT + torS -> OMG clustered90 reps (~660M reps,
# 171 GB, ~10x the 10M pilot). mmseqs unsplit memory footprint would be
# ~1.7 TB, so --split-memory-limit 200G forces sequential DB splits that fit
# inside the 256G request.
#
# --max-seqs 200000 replaces the pilot's 5000 cap. The pilot saturated at
# 4999/5000 torS hits; 200k gives comfortable headroom even scaled 10x.

set -euo pipefail

REPO=/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/feder_harris_hackathon_26
DB=/fh/working/srivatsan_s/databases/mmseqs2/omg_cds/omg_cds_clustered90_seqs

cd "$REPO"
mkdir -p results tmp logs

ml MMseqs2/13-45111-gompi-2021b
mmseqs version

mmseqs createdb queries/torT_torS.faa queries/torT_torS_db

mmseqs search \
    queries/torT_torS_db \
    "$DB" \
    results/torT_torS_vs_omgfull \
    tmp/torT_torS_full \
    -s 7.5 \
    -e 1e-10 \
    --max-seqs 200000 \
    --split-memory-limit 200G \
    --threads "$SLURM_CPUS_PER_TASK"

mmseqs convertalis \
    queries/torT_torS_db \
    "$DB" \
    results/torT_torS_vs_omgfull \
    results/torT_torS_vs_omgfull.m8 \
    --format-output "query,target,theader,pident,alnlen,evalue,bits,qcov,tcov"

echo "done. hits: $(wc -l < results/torT_torS_vs_omgfull.m8)"
