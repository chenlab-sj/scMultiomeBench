#!/bin/bash

#BSUB -P benchmark
#BSUB -J fig7_sumplot
#BSUB -q standard
#BSUB -R rusage[mem=16000]
#BSUB -n 1
#BSUB -o sum_perfom_plot.log
#BSUB -e sum_perfom_plot.err

# Fig7 cross-dataset summary: regenerate sum_perform.pdf + sum_perform_rank.pdf.
# The only figure change is the rank-plot method-color legend, now forced to 2 rows (guide_legend(nrow = 2)).
# sum_perfom_plot.R auto-detects cluster vs Mac paths, so on the cluster it reads inputs from ${TOOLS_ROOT}.
#
# R packages required: dplyr, tidyr, purrr, RColorBrewer, ggplot2 (all in seurat4).
# ggbump is NO LONGER needed -- the rank-plot bump curves are now drawn by a self-contained
# sigmoid helper inside sum_perfom_plot.R, so nothing extra has to be installed.
module load conda3/202311
conda activate seurat4
cd ${PROJECT_ROOT}/sum_plot

Rscript sum_perfom_plot.R
RC=$?
conda deactivate

if [ "$RC" -eq 0 ]; then
  echo "DONE: sum_perform.pdf + sum_perform_rank.pdf (rank legend = 2 rows)"
else
  echo "FAILED: Rscript exited with code $RC -- see sum_perfom_plot.err (likely a missing R package)"
fi
