#!/bin/bash

#BSUB -P benchmark
#BSUB -J seurat3
#BSUB -q compbio
#BSUB -R rusage[mem=100000]
#BSUB -n 1
#BSUB -o seurat3.log
#BSUB -e seurat3.err


module load conda3/202210 

conda activate seurat4

R CMD BATCH seurat3.R

conda deactivate

