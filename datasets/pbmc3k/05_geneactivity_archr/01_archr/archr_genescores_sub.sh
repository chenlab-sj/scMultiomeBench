#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P archr
#BSUB -J pbmc3k_archr_genescores
#BSUB -q large_mem
#BSUB -R rusage[mem=200000]
#BSUB -n 8
#BSUB -o archr_genescores.log
#BSUB -e archr_genescores.err

# ArchR gene scores for pbmc3k (R2.1 alt gene-activity). seurat4 env has ArchR.
module load conda3/202210
conda activate seurat4

cd ${PROJECT_ROOT}/pbmc/pbmc3k/benchmark/geneactivity_archr/archr

R CMD BATCH archr_genescores_pbmc3k.R

conda deactivate
echo "DONE: ArchR_pbmc3k/export/archr_gene_scores.mtx (+ genes.csv, barcodes.csv)"
