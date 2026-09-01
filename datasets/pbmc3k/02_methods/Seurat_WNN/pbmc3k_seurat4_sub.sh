#!/bin/bash

#BSUB -P seurat4
#BSUB -J pbmc3k_seurat4
#BSUB -q standard
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o pbmc3k_seurat4.log
#BSUB -e pbmc3k_seurat4.err


module load conda3/202210 
conda activate seurat4

R CMD BATCH pbmc3k_seurat4.R

conda deactivate
