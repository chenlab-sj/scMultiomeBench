#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J parse_plot
#BSUB -q large_mem
#BSUB -R rusage[mem=40000]
#BSUB -n 1
#BSUB -o parse_plot.log
#BSUB -e parse_plot.err

# Builds metrics_matrix.csv + sum_metrics_clean.csv + metrics_table.html for the Parse cross-platform.
# SPLICE_PUBLISHED=0 -> score every method from THIS run (standalone Parse figure; no pbmc3k merge).
# Needs R with dplyr + formattable (+ webshot2/webshot for the PNG). Can also just run on your Mac:
#     SPLICE_PUBLISHED=0 Rscript plot_metrics_matrix.R . .
module load conda3/202311
conda activate seurat4          # <-- adjust to an R env that has dplyr + formattable

ROOT=${PROJECT_ROOT}
cd "${ROOT}/pbmc_parse/benchmark/metrics"

export SPLICE_PUBLISHED=0
Rscript plot_metrics_matrix.R . .

conda deactivate
echo "DONE: metrics_matrix.csv + sum_metrics_clean.csv + metrics_table.html"
