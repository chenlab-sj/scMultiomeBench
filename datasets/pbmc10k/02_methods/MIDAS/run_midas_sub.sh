#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J run_midas
#BSUB -q gpu_interactive
#BSUB -R rusage[mem=100000]
#BSUB -gpu "num=1"
#BSUB -o run_midas.log
#BSUB -e run_midas.err

module load conda3/202311

source activate midas-env

export PYTHONNOUSERSITE=1

cd ${PROJECT_ROOT}/pbmc/pbmc10k/scripts/midas

python 01_run_midas.py

conda deactivate
