#!/bin/bash

#BSUB -P benchmark
#BSUB -J scBridge_s0
#BSUB -q standard
#BSUB -R rusage[mem=100000]
#BSUB -n 2
#BSUB -o scBridge_s0.log
#BSUB -e scBridge_s0.err


module load conda3/202210 

conda activate seurat4

R CMD BATCH 00_scBridge_s0.R

conda deactivate

