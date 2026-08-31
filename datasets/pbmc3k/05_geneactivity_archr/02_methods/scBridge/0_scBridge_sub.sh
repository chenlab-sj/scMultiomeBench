#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P scBridge
#BSUB -J scBridge_s0_archr
#BSUB -q standard
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o scBridge_s0_archr.log
#BSUB -e scBridge_s0_archr.err


module load conda3/202210
conda activate seurat4

cd ${PROJECT_ROOT}/pbmc/pbmc3k/benchmark/geneactivity_archr/scBridge/
R CMD BATCH 00_make_scbridge_h5ad.R
