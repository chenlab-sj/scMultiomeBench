#!/bin/bash

#BSUB -P benchmark
#BSUB -J seurat3
#BSUB -q standard
#BSUB -R rusage[mem=200000]
#BSUB -n 1
#BSUB -o seurat3.log
#BSUB -e seurat3.err


module load conda3/202210 

conda activate seurat4

R CMD BATCH macro_subset_rep1.R

conda deactivate

