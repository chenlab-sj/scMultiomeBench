#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J rms_run_maxfuse
#BSUB -q large_mem
#BSUB -R rusage[mem=100000]
#BSUB -n 4
#BSUB -o run_maxfuse.log
#BSUB -e run_maxfuse.err

# MaxFuse Fusor rep1 (CPU). Run AFTER 00_prep_maxfuse_input_sub.sh.
module load conda3/202311
source activate maxfuse-env
export PYTHONNOUSERSITE=1
cd ${PROJECT_ROOT}/RMS/Mast607/script/maxfuse
python 01_run_maxfuse.py
conda deactivate
echo "DONE: latent.csv"
