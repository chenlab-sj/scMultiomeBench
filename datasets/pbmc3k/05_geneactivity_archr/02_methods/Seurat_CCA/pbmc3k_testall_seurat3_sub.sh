#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P Seurat3
#BSUB -J pbmc3k_testall_seurat3_archr
#BSUB -q compbio
#BSUB -R rusage[mem=100000]
#BSUB -n 2
#BSUB -o pbmc3k_testall_seurat3_archr.log
#BSUB -e pbmc3k_testall_seurat3_archr.err


module load conda3/202210

conda activate seurat4

cd ${PROJECT_ROOT}/pbmc/pbmc3k/benchmark/geneactivity_archr/Seurat_CCA
mkdir -p res_pbmc3k_testall

R CMD BATCH pbmc3k_testall_seurat3.R

conda deactivate
