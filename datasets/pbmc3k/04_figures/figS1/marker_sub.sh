#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J pbmc3k_marker_plot
#BSUB -q standard
#BSUB -R rusage[mem=16000]
#BSUB -n 1
#BSUB -o marker_plot.log
#BSUB -e marker_plot.err

# R1.6 label-validity marker plots: canonical markers x cell type, faceted by the cell type each gene marks.
# Rows = the GROUND-TRUTH labels used to score the benchmark (the labels under scrutiny); the SingleR run is
# an optional cross-check. NB plot_marker_label_validation.R auto-detects cluster vs Mac paths and runs fine locally
# (just `Rscript plot_marker_label_validation.R`) -- this job is only needed if Seurat isn't available on your machine.
module load conda3/202311
conda activate seurat4
cd ${PROJECT_ROOT}/pbmc/pbmc3k/benchmark/label_sensitivity
Rscript plot_marker_label_validation.R                        # 10X ground-truth labels -> marker_plot_label.pdf
LABELS=singleR Rscript plot_marker_label_validation.R         # SingleR RNA-only labels  -> marker_plot_singleR.pdf
conda deactivate
echo "DONE: marker_plot_label.pdf + marker_plot_singleR.pdf"
