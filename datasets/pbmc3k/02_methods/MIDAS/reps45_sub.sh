#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P benchmark
#BSUB -J pbmc3k_midas_reps45
#BSUB -q gpu_interactive
#BSUB -R rusage[mem=100000]
#BSUB -gpu "num=1"
#BSUB -o run_midas_reps45.log
#BSUB -e run_midas_reps45.err
# Two ADDITIONAL reproducibility reps for MIDAS (seeds 100, 200; rep1-3 used 420/0/42) ->
# midas/rep4/latent.csv, midas/rep5/latent.csv. 01_run_midas.py reads REP+SEED from env and writes to
# scripts/midas/$REP/. pbmc3k is small, so the two trainings run serially in one job.
module load conda3/202311
source activate midas-env
export PYTHONNOUSERSITE=1
cd ${PROJECT_ROOT}/pbmc/pbmc3k/scripts/midas
for rs in "rep4 100" "rep5 200"; do
  set -- $rs; export REP=$1 SEED=$2
  echo "=== MIDAS $REP (seed $SEED) ==="
  python 01_run_midas.py
done
conda deactivate
echo "DONE: midas/rep4/latent.csv + midas/rep5/latent.csv"
