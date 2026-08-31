#!/bin/bash

#BSUB -P benchmark
#BSUB -J bmmc_sp_plot
#BSUB -q large_mem
#BSUB -R rusage[mem=40000]
#BSUB -n 1
#BSUB -o cs_plot.log
#BSUB -e cs_plot.err

# BMMC cross-site / paired matrix from the metric CSVs. SPLICE_PUBLISHED=0 -> score+rank ALL methods
# from THIS run (major/ basis) and DROP the same-cell columns (ks.statistic/ari/ami) so the crosssite
# (unpaired) and s1d1_paired (multiome) matrices are scored on the identical 6 metrics -> fair scatter.
# Runs after major-recompute AND peak. LS_SUBCWD = the experiment dir.
module load conda3/202311
conda activate seurat4
cd "${LS_SUBCWD:-$PWD}" || exit 1
export SPLICE_PUBLISHED=0
Rscript plot_metrics_matrix.R . .
conda deactivate
echo "DONE: metrics_matrix.csv + sum_metrics_clean.csv + metrics_table.html/png"
