#!/bin/bash

#BSUB -P Seurat3
#BSUB -J pbmc3k_testall_seurat3
#BSUB -q rhel8_compbio
#BSUB -R rusage[mem=100000]
#BSUB -n 2
#BSUB -o pbmc3k_testall_seurat3.log
#BSUB -e pbmc3k_testall_seurat3.err


module load conda3/202210 

conda activate seurat4

R CMD BATCH pbmc3k_testall_seurat3.R

conda deactivate
