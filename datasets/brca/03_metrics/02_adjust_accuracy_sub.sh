#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J brca_figS4a_adjaccu
#BSUB -q standard
#BSUB -R rusage[mem=20000]
#BSUB -n 1
#BSUB -o figS4a_adjaccu.log
#BSUB -e figS4a_adjaccu.err

# Random-adjusted ATAC cell-type prediction accuracy (average_accu) -> adj_atac_predaccu.csv.
# Reads all knn_pred_label__*.csv in this dir, so it reflects whatever methods are currently
# present (14 after the Conos/LIGER cleanup). Fast; benchmark_env (sklearn).
module load conda3/202303
conda activate benchmark_env

ROOT=${PROJECT_ROOT}
cd "${ROOT}/BRCA/benchmark/figS4a"

python 01_adjust_accuracy.py --label label.csv --out .

conda deactivate
echo "DONE: adj_atac_predaccu.csv"
