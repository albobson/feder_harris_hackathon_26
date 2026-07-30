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
#$ -l h_rt=2:0:0
#$ -cwd
#$ -j y
#$ -o jeffrey/logs/ortholog_map.$JOB_ID.log

set -euo pipefail

echo "host      : $(hostname)"
echo "job id    : ${JOB_ID:-none}"
echo "slots     : ${NSLOTS:-1}"
echo "started   : $(date)"
echo

# NOTE: do not use $(dirname "$0") here -- under SGE, $0 is the spooled copy of
# this script in /var/spool/uge/..., not its path in the repo. Because -cwd is
# set, the job starts in the directory qsub was run from (the repo root), which
# SGE also exposes as $SGE_O_WORKDIR.
cd "${SGE_O_WORKDIR:-$PWD}/jeffrey/analysis"
python3 build_ortholog_map.py

echo
echo "finished  : $(date)"
