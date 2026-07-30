#!/bin/bash
# Submit the BLAST-based ortholog mapping to an SGE compute node.
# blastp across 9 proteomes is ~CPU-hours, so it must not run on a login node.
#
#   qsub jeffrey/analysis/submit_ortholog_map.sh
#
#$ -N ortholog_map
#$ -q feder-short.q
#$ -pe serial 8
#$ -l mfree=4G
#$ -cwd
#$ -j y
#$ -o jeffrey/logs/ortholog_map.$JOB_ID.log

set -euo pipefail

echo "host      : $(hostname)"
echo "job id    : ${JOB_ID:-none}"
echo "slots     : ${NSLOTS:-1}"
echo "started   : $(date)"
echo

cd "$(dirname "$0")"
python3 build_ortholog_map.py

echo
echo "finished  : $(date)"
