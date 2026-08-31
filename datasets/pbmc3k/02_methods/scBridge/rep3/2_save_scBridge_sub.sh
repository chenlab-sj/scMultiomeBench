#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P scBridge
#BSUB -J scBridge_s2
#BSUB -q standard
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o scBridge_s2.log
#BSUB -e scBridge_s2.err


module load conda3/202210

conda activate scBridge
module load gcc/13.1.0-rhel7

cd ${TOOLS_ROOT}/scBridge
python 2_save_scBridge_result.py

conda deactivate
