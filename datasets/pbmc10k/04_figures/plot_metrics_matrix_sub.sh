#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J figS1c_plot
#BSUB -q large_mem
#BSUB -R rusage[mem=40000]
#BSUB -n 1
#BSUB -o figS1c_plot.log
#BSUB -e figS1c_plot.err

# FigS1C composite table: build fig2b_matrix.csv + sum_metrics_clean.csv + fig2b_table.html (+ PNG/PDF if the
# env has webshot2/Chrome) from the 23-method metric CSVs (metrics + peak jobs must be done first).
# SPLICE_PUBLISHED=0 = score ALL 23 methods from this run (full recompute; peak already spliced old+new in
# peakdist.csv). seurat4 = R with dplyr + formattable. If the PNG/PDF doesn't render here (no webshot2),
# the HTML is still written -> render it to PNG on the Mac via headless Chrome (same as figS4a).
module load conda3/202311
conda activate seurat4

ROOT=${PROJECT_ROOT}
cd "${ROOT}/pbmc/pbmc10k/benchmark"

export SPLICE_PUBLISHED=0
Rscript plot_metrics_matrix.R . .

conda deactivate
echo "DONE: fig2b_matrix.csv + sum_metrics_clean.csv + fig2b_table.html"
