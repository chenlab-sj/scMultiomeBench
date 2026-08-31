#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J parse_adjaccu
#BSUB -q large_mem
#BSUB -R rusage[mem=16000]
#BSUB -n 1
#BSUB -o parse_adjaccu.log
#BSUB -e parse_adjaccu.err

# Random-adjusted ATAC cell-type prediction accuracy (average_accu) -> adj_atac_predaccu.csv.
# Combines ALL knn_pred_label__*.csv written by 00_compute_metrics_sub.sh with label.csv's random_atac
# baseline. RUN AFTER 00_compute_metrics_sub.sh finishes (needs every method's knn_pred_label__<m>.csv).
# Light + idempotent -> also fine to re-run standalone if you add a method later.
module load conda3/202303
conda activate benchmark_env

ROOT=${PROJECT_ROOT}
cd "${ROOT}/pbmc_parse/benchmark/metrics"

python 01_adjust_accuracy.py --label label.csv --out .

conda deactivate
echo "DONE: adj_atac_predaccu.csv"
