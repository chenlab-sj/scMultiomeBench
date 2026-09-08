#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P benchmark
#BSUB -J ht137_scbutterfly_rep2
#BSUB -q dgx
#BSUB -R rusage[mem=200000]
#BSUB -gpu "num=1"
#BSUB -o run_scbutterfly_rep2.log
#BSUB -e run_scbutterfly_rep2.err

# HT137B1-S1H7 scButterfly rep2 (seed 100), FORCE_SEED=1 to override the internal 19193 -> genuine seed
# variation. Single training. base latent.csv untouched (REP="" only).
module load conda3/202311
source activate scbutterfly-env
export PYTHONNOUSERSITE=1
export FORCE_SEED=1

cd ${PROJECT_ROOT}/BRCA/HT137B1-S1H7/scbutterfly
export REP=rep2 SEED=100
echo "=== scButterfly rep2 (seed 100, FORCE_SEED=1) ==="
python run_scbutterfly.py

conda deactivate
echo "DONE: HT137B1-S1H7/scbutterfly/rep2/latent.csv"
