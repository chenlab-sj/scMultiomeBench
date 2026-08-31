#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P benchmark
#BSUB -J azimuth_parse
#BSUB -q standard            # <-- EDIT to a CPU queue on your cluster (you used large_mem for BABEL).
#BSUB -n 4
#BSUB -R rusage[mem=100000]
#BSUB -o azimuth_parse.log
#BSUB -e azimuth_parse.err

# Annotate ONE Parse donor with Azimuth (reproduces Parse Fig 1), harmonize to the
# 6-type taxonomy, export label CSV + UMAPs. Change DONOR and resubmit for each donor.
#
# !!! INTERNET: RunAzimuth downloads the "pbmcref" reference on first use. Compute nodes
#     are often offline (the same trap that broke BABEL's leidenalg/igraph download).
#     EITHER run this on a node WITH internet, OR pre-stage pbmcref once on a login node:
#         module load conda3/202311; source activate seurat4
#         Rscript -e 'library(SeuratData); InstallData("pbmcref")'
#     After it's cached, this batch job can reuse it offline.

DONOR="Donor_4"

module load conda3/202210 

source activate seurat4     # <-- EDIT: env with Seurat + Azimuth + SeuratData(pbmcref)
export PYTHONNOUSERSITE=1

cd ${PROJECT_ROOT}/pbmc_parse/data
Rscript 01_azimuth_annotate_parse.R "$DONOR"

conda deactivate
