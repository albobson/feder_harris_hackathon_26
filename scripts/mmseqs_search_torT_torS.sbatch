#!/bin/bash
#SBATCH --job-name=mmseqs_torT_torS
#SBATCH --partition=campus-new
#SBATCH --gpus=rtx2080ti:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/mmseqs_torT_torS_%j.out
#SBATCH --error=logs/mmseqs_torT_torS_%j.err

# GPU mmseqs2 search: MG1655 torT + torS -> OMG 10M padded cluster DB.
# 18 GB padded DB fits comfortably on an RTX 2080ti (11 GB VRAM) because
# mmseqs streams k-mer prefiltering; ungapped/gapped alignment blocks are
# small. If VRAM turns out to be the bottleneck, resubmit against chorus:
#   #SBATCH --partition=chorus
#   #SBATCH --gpus=l40s:1
# (note: chorus uses a separate module tree — may need to build MMseqs2
#  from source or use apptainer/singularity there.)

set -euo pipefail

REPO=/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/feder_harris_hackathon_26
DB=/fh/working/srivatsan_s/databases/mmseqs2/omg_cds_10M/omg_cds_clustered90_seqs_pad

cd "$REPO"
mkdir -p results tmp logs

ml MMseqs2/13-45111-gompi-2021b
mmseqs version
nvidia-smi

mmseqs createdb queries/torT_torS.faa queries/torT_torS_db

mmseqs search \
    queries/torT_torS_db \
    "$DB" \
    results/torT_torS_vs_omg10M \
    tmp/torT_torS \
    --gpu 1 \
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
