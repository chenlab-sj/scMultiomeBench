#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J brca_cobolt_nmi
#BSUB -q large_mem
#BSUB -R rusage[mem=100000]
#BSUB -n 1
#BSUB -o cobolt_nmi.log
#BSUB -e cobolt_nmi.err

# Cluster all 5 BRCA Cobolt reps and report 3-rep (current figure) vs 5-rep (all seeds) NMI -> cobolt_nmi.log.
module load conda3/202303
conda activate benchmark_env
ROOT=${PROJECT_ROOT}
cd "${ROOT}/BRCA/benchmark/fig4"
export BENCHMARK_FUN_DIR=${BENCHMARK_FUN_DIR}
export LABEL=label.csv
export FIX_RES=0
python 00_cobolt_nmi.py
conda deactivate
