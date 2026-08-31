#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P benchmark
#BSUB -J ga_archr_plot
#BSUB -q large_mem
#BSUB -R rusage[mem=40000]
#BSUB -n 1
#BSUB -o ga_archr_plot.log
#BSUB -e ga_archr_plot.err
# PAIRED composite (pbmc3k same cells) -> plot_metrics_matrix.R, full metrics. Score only the 6 ArchR variants.
module load conda3/202311
conda activate seurat4
cd ${PROJECT_ROOT}/pbmc/pbmc3k/benchmark/geneactivity_archr/metrics
export SPLICE_PUBLISHED=0
Rscript plot_metrics_matrix.R . .
conda deactivate
echo "DONE: fig2b_matrix.csv + sum_metrics_clean.csv (ArchR variant scores)"
