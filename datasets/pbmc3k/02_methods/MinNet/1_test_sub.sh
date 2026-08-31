#!/bin/bash

#BSUB -P MinNet
#BSUB -J MinNet_test1
#BSUB -q rhel8_compbio
#BSUB -R rusage[mem=100000]
#BSUB -n 2
#BSUB -o MinNet_test1.log
#BSUB -e MinNet_test1.err


module load conda3/202210 

conda activate seurat4

R CMD BATCH 1_test_dataprocess.R

conda deactivate   
