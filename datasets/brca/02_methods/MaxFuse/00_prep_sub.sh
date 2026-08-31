#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J brca_maxfuse_prep
#BSUB -q large_mem
#BSUB -R rusage[mem=200000]
#BSUB -n 1
#BSUB -o prep_maxfuse.log
#BSUB -e prep_maxfuse.err

# Signac gene-activity prep for BRCA. seurat4 env (Seurat/Signac + reticulate->anndata).
# large_mem + 200G: the BRCA ATAC assay (14G h5ad) + GeneActivity is far heavier than pbmc3k.
module load conda3/202210
conda activate seurat4

cd ${PROJECT_ROOT}/BRCA/HT243B1-S1H4/maxfuse

R CMD BATCH 00_prep_maxfuse_input.R

conda deactivate
echo "DONE: 4 input h5ads in maxfuse/input/"
