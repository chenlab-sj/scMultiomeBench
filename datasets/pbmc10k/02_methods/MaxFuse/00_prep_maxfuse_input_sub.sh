#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J maxfuse_prep
#BSUB -q compbio
#BSUB -R rusage[mem=100000]
#BSUB -n 1
#BSUB -o prep_maxfuse.log
#BSUB -e prep_maxfuse.err

# Signac gene-activity prep. Use the `seurat4` env (Seurat/Signac + reticulate->anndata
# for write_h5ad), exactly like 00_scBridge_s0.R / 00_scjoint_s1.R.
module load conda3/202210

conda activate seurat4

cd ${PROJECT_ROOT}/pbmc/pbmc10k/scripts/maxfuse

R CMD BATCH 00_prep_maxfuse_input.R

conda deactivate
