#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P benchmark
#BSUB -J brca_run_midas
#BSUB -q dgx
#BSUB -R rusage[mem=200000]
#BSUB -gpu "num=1"
#BSUB -o run_midas.log
#BSUB -e run_midas.err

# MIDAS (GPU). 200G mem: the BRCA ATAC (14G h5ad) + train concat is far heavier than pbmc3k.
# Smoke first if you like: export MIDAS_EPOCHS=10 before bsub to check the pipeline end-to-end.
module load conda3/202311
source activate midas-env
export PYTHONNOUSERSITE=1

cd ${PROJECT_ROOT}/BRCA/HT137B1-S1H7/midas

python run_midas.py

conda deactivate
echo "DONE: latent.csv"
