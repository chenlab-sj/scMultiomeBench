#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J bmmc_run_midas
#BSUB -q dgx
#BSUB -R rusage[mem=150000]
#BSUB -gpu "num=1"
#BSUB -o run_midas.log
#BSUB -e run_midas.err

# MIDAS (WITHOUT batch) on BMMC_d1 (GPU). Smoke: export MIDAS_EPOCHS=10 before bsub.
module load conda3/202311
source activate midas-env
export PYTHONNOUSERSITE=1
cd ${PROJECT_ROOT}/BMMC_d1/scripts/crosssite/midas
python 01_run_midas.py
conda deactivate
echo "DONE: latent.csv"
