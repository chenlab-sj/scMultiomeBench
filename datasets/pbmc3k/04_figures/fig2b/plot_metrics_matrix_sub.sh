#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J fig2b_plot
#BSUB -q large_mem
#BSUB -R rusage[mem=40000]
#BSUB -n 1
#BSUB -o fig2b_plot.log
#BSUB -e fig2b_plot.err

# Builds fig2b_matrix.csv + sum_metrics_clean.csv + fig2b_table.html from the metric CSVs.
# Needs R with dplyr + formattable (+ htmlwidgets/webshot for the PNG). Swap the env below
# to whichever of your R envs has them (seurat4 / bindsc), or just run on your Mac:
#     Rscript plot_metrics_matrix.R . .
module load conda3/202311
conda activate seurat4          # <-- adjust to an R env that has dplyr + formattable

ROOT=${PROJECT_ROOT}
cd "${ROOT}/pbmc/pbmc3k/benchmark/fig2b"

Rscript plot_metrics_matrix.R . .

conda deactivate
echo "DONE: fig2b_matrix.csv + fig2b_table.html"
