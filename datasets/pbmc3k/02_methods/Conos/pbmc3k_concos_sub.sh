#!/bin/bash

#BSUB -P concos
#BSUB -J pbmc3k_concos
#BSUB -q rhel8_standard
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o pbmc3k_concos.log
#BSUB -e pbmc3k_concos.err


module load conda3/202210 
conda activate seurat4

R CMD BATCH pbmc3k_concos.R

conda deactivate
