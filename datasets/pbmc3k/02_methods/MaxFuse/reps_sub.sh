#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J pbmc3k_maxfuse_reps
#BSUB -q large_mem
#BSUB -R rusage[mem=100000]
#BSUB -n 4
#BSUB -o run_maxfuse_reps.log
#BSUB -e run_maxfuse_reps.err

# MaxFuse rep2 + rep3 (seeds 0, 42) -> maxfuse/rep2/latent.csv, maxfuse/rep3/latent.csv.
# Reuses the shared maxfuse/input/ from the rep1 prep (deterministic) -> NO prep re-run.
module load conda3/202311
source activate maxfuse-env
export PYTHONNOUSERSITE=1

cd ${PROJECT_ROOT}/pbmc/pbmc3k/scripts/maxfuse

for rs in "rep2 0" "rep3 42"; do
  set -- $rs; export REP=$1 SEED=$2
  echo "=== MaxFuse $REP (seed $SEED) ==="
  python 01_run_maxfuse.py
done

conda deactivate
echo "DONE: maxfuse/rep2/latent.csv + maxfuse/rep3/latent.csv"
