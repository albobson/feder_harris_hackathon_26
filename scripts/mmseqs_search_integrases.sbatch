#!/bin/bash
#SBATCH --job-name=mmseqs_integrases
#SBATCH --partition=campus-new
#SBATCH --gpus=rtx2080ti:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=logs/mmseqs_integrases_%j.out
#SBATCH --error=logs/mmseqs_integrases_%j.err

# GPU mmseqs2 search: 11 phage integrases -> OMG 10M padded cluster DB.
# See mmseqs_search_torT_torS.sbatch for chorus/L40S fallback details.

set -euo pipefail

REPO=/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/feder_harris_hackathon_26
DB=/fh/working/srivatsan_s/databases/mmseqs2/omg_cds_10M/omg_cds_clustered90_seqs_pad

cd "$REPO"
mkdir -p results tmp logs

ml MMseqs2/13-45111-gompi-2021b
mmseqs version
nvidia-smi

mmseqs createdb queries/phage_integrases.faa queries/phage_integrases_db

mmseqs search \
    queries/phage_integrases_db \
    "$DB" \
    results/integrases_vs_omg10M \
    tmp/integrases \
    --gpu 1 \
    -s 7.5 \
    -e 1e-10 \
    --max-seqs 20000 \
    --threads "$SLURM_CPUS_PER_TASK"

mmseqs convertalis \
    queries/phage_integrases_db \
    "$DB" \
    results/integrases_vs_omg10M \
    results/integrases_vs_omg10M.m8 \
    --format-output "query,target,theader,pident,alnlen,evalue,bits,qcov,tcov"

echo "done. hits: $(wc -l < results/integrases_vs_omg10M.m8)"
