#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J pbmc3k_singleR
#BSUB -q standard
#BSUB -R rusage[mem=32000]
#BSUB -n 1
#BSUB -o singleR.log
#BSUB -e singleR.err

# Independent RNA-only re-annotation of pbmc3k with SingleR (R4/R1.6 label-circularity rebuttal).
# PREREQ: the MonacoImmuneData reference must be reachable. If compute nodes have no internet, pre-fetch
# once on a login node:  module load conda3/202311; conda activate seurat4; R -e 'celldex::MonacoImmuneData()'
module load conda3/202311
conda activate seurat4
cd ${PROJECT_ROOT}/pbmc/pbmc3k/benchmark/label_sensitivity
Rscript 00_singleR_annotate.R
conda deactivate
echo "DONE: label_singleR.csv + confusion_10x_vs_singleR.csv + singleR_per_cell.csv"
