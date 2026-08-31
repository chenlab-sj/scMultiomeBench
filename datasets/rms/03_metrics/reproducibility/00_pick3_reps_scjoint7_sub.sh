#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J rms_pick3_scjoint7
#BSUB -q large_mem
#BSUB -R rusage[mem=100000]
#BSUB -n 1
#BSUB -o pick3_scjoint7.log
#BSUB -e pick3_scjoint7.err

# Run after RMS/Mast607/script/scJoint/reps67_sub.sh finishes. This restages
# scJoint reps 1-7 and evaluates all 3-rep subsets into pick3_scjoint7/.
set -euo pipefail

module load conda3/202303
conda activate benchmark_env

ROOT=${PROJECT_ROOT}

cd "${ROOT}/RMS/benchmark"
python 01_prep_latents.py

cd "${ROOT}/RMS/benchmark/fig4"
export BENCHMARK_FUN_DIR=${BENCHMARK_FUN_DIR}
export LABEL=label.csv
export FIX_RES="${FIX_RES:-0}"
export METHODS=scJoint
export OUT=pick3_scjoint7

python 00_pick3_reps.py

conda deactivate
