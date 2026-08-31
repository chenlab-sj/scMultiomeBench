#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J brca_midas_pick
#BSUB -q large_mem
#BSUB -R rusage[mem=100000]
#BSUB -n 1
#BSUB -o midas_pick.log
#BSUB -e midas_pick.err

# Cluster all 5 BRCA MIDAS reps and report Fig4A repro + Fig4B NMI for each 3-rep subset (-> midas_pick.csv).
# FIX_RES=0 matches the current figure; set FIX_RES=1 to also see the fixed-resolution NMI.
module load conda3/202303
conda activate benchmark_env

ROOT=${PROJECT_ROOT}
cd "${ROOT}/BRCA/benchmark/fig4"
export BENCHMARK_FUN_DIR=${BENCHMARK_FUN_DIR}
export LABEL=label.csv
export FIX_RES=0

python 01_midas_pick.py

conda deactivate
echo "DONE: midas_pick.csv (each MIDAS 3-subset's Fig4A repro + Fig4B NMI vs Cobolt 0.6746)"
