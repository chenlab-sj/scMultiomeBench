#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P benchmark
#BSUB -J brca_run_maxfuse
#BSUB -q large_mem
#BSUB -R rusage[mem=100000]
#BSUB -n 4
#BSUB -o run_maxfuse.log
#BSUB -e run_maxfuse.err

# MaxFuse Fusor (CPU-only). Reads the 4 small input h5ads from prep -> latent.csv.
# Run AFTER prep_sub.sh finishes.
module load conda3/202311
source activate maxfuse-env

# Isolate from ~/.local user-site so the env's packages aren't shadowed.
export PYTHONNOUSERSITE=1

cd ${PROJECT_ROOT}/BRCA/HT137B1-S1H7/maxfuse

python run_maxfuse.py

conda deactivate
echo "DONE: latent.csv"
