#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J run_mira_unp
#BSUB -q dgx
#BSUB -R rusage[mem=100000]
#BSUB -gpu "num=1"
#BSUB -o run_mira_unpaired.log
#BSUB -e run_mira_unpaired.err

module load conda3/202311

source activate mira-env

# Isolate from ~/.local user-site so the env's mira-multiome isn't shadowed.
export PYTHONNOUSERSITE=1

cd ${PROJECT_ROOT}/pbmc/pbmc10k/scripts/MIRA

python run_mira_unpaired.py

conda deactivate
