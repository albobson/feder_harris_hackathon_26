#!/bin/bash
#SBATCH --job-name=mmseqs_integrases
#SBATCH --partition=short
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=logs/mmseqs_integrases_%j.out
#SBATCH --error=logs/mmseqs_integrases_%j.err

# CPU mmseqs2 search: 10 phage integrases -> OMG 10M cluster reps.
# See mmseqs_search_torT_torS.sbatch for the v13/GPU note.

set -euo pipefail

REPO=/fh/fast/srivatsan_s/grp/SrivatsanLab/Dustin/feder_harris_hackathon_26
DB=/fh/working/srivatsan_s/databases/mmseqs2/omg_cds_10M/omg_cds_clustered90_seqs

cd "$REPO"
mkdir -p results tmp logs

ml MMseqs2/13-45111-gompi-2021b
mmseqs version

mmseqs createdb queries/phage_integrases.faa queries/phage_integrases_db

mmseqs search \
    queries/phage_integrases_db \
    "$DB" \
    results/integrases_vs_omg10M \
    tmp/integrases \
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
