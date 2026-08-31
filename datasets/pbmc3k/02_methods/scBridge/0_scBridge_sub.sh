#!/bin/bash

#BSUB -P scBridge
#BSUB -J scBridge_s0
#BSUB -q rhel8_standard
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o scBridge_s0.log
#BSUB -e scBridge_s0.err


module load conda3/202210 
conda activate seurat4

R CMD BATCH 00_make_scbridge_h5ad.R
