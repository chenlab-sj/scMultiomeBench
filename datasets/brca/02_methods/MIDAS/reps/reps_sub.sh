#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J brca_midas_reps
#BSUB -q dgx
#BSUB -R rusage[mem=200000]
#BSUB -gpu "num=1"
#BSUB -o run_midas_reps.log
#BSUB -e run_midas_reps.err

# MIDAS rep2 + rep3 (seeds 0, 42) -> midas/rep2/latent.csv, midas/rep3/latent.csv.
# Two full GPU trainings run SERIALLY here. If the queue time limit is tight, split into two
# jobs: set REP/SEED to one pair, drop the loop, and submit twice (they run in parallel).
module load conda3/202311
source activate midas-env
export PYTHONNOUSERSITE=1

cd ${PROJECT_ROOT}/BRCA/HT243B1-S1H4/midas

for rs in "rep2 0" "rep3 42"; do
  set -- $rs; export REP=$1 SEED=$2
  echo "=== MIDAS $REP (seed $SEED) ==="
  python 01_run_midas.py
done

conda deactivate
echo "DONE: midas/rep2/latent.csv + midas/rep3/latent.csv"
