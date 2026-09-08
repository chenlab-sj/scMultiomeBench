#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P benchmark
#BSUB -J ht137_scbutterfly_rep3
#BSUB -q dgx
#BSUB -R rusage[mem=200000]
#BSUB -gpu "num=1"
#BSUB -o run_scbutterfly_rep3.log
#BSUB -e run_scbutterfly_rep3.err

# HT137B1-S1H7 scButterfly rep3 (seed 200), FORCE_SEED=1. Single training. base + rep2 untouched.
module load conda3/202311
source activate scbutterfly-env
export PYTHONNOUSERSITE=1
export FORCE_SEED=1

cd ${PROJECT_ROOT}/BRCA/HT137B1-S1H7/scbutterfly
export REP=rep3 SEED=200
echo "=== scButterfly rep3 (seed 200, FORCE_SEED=1) ==="
python run_scbutterfly.py

conda deactivate
echo "DONE: HT137B1-S1H7/scbutterfly/rep3/latent.csv"
