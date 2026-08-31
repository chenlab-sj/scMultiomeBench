#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J rms_pick3_reps
#BSUB -q large_mem
#BSUB -R rusage[mem=100000]
#BSUB -n 1
#BSUB -o pick3_reps.log
#BSUB -e pick3_reps.err

# Evaluate every 3-rep subset for selected RMS Fig4 methods.
# Override METHODS if needed, e.g.:
#   METHODS=scJoint,scVI,Portal,BindSC,Cobolt,scDART bsub < 00_pick3_reps_sub.sh
module load conda3/202303
conda activate benchmark_env

ROOT=${PROJECT_ROOT}
cd "${ROOT}/RMS/benchmark/fig4"
export BENCHMARK_FUN_DIR=${BENCHMARK_FUN_DIR}
export LABEL=label.csv
export FIX_RES="${FIX_RES:-0}"
export METHODS="${METHODS:-scJoint,scVI,Portal,BindSC,Cobolt,scDART}"

python 00_pick3_reps.py

conda deactivate
