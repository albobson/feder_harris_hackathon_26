#!/bin/bash
#SBATCH --job-name=mmseqs_torT_torS
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/mmseqs_torT_torS_%j.out
#SBATCH --error=logs/mmseqs_torT_torS_%j.err

# CPU mmseqs2 search: MG1655 torT + torS -> OMG 10M cluster reps.
# GPU search would need MMseqs2 v14+; the cluster only has v13, so we
# fall back to CPU on the non-padded cluster reps.

set -euo pipefail

REPO=/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/feder_harris_hackathon_26
DB=/fh/working/srivatsan_s/databases/mmseqs2/omg_cds_10M/omg_cds_clustered90_seqs

cd "$REPO"
mkdir -p results tmp logs

ml MMseqs2/13-45111-gompi-2021b
mmseqs version

mmseqs createdb queries/torT_torS.faa queries/torT_torS_db

mmseqs search \
    queries/torT_torS_db \
    "$DB" \
    results/torT_torS_vs_omg10M \
    tmp/torT_torS \
    -s 7.5 \
    -e 1e-10 \
    --max-seqs 5000 \
    --threads "$SLURM_CPUS_PER_TASK"

mmseqs convertalis \
    queries/torT_torS_db \
    "$DB" \
    results/torT_torS_vs_omg10M \
    results/torT_torS_vs_omg10M.m8 \
    --format-output "query,target,theader,pident,alnlen,evalue,bits,qcov,tcov"

echo "done. hits: $(wc -l < results/torT_torS_vs_omg10M.m8)"
