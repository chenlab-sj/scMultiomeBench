#!/bin/bash

#BSUB -P benchmark
#BSUB -J scjoint_s1
#BSUB -q rhel8_compbio
#BSUB -R rusage[mem=100000]
#BSUB -n 2
#BSUB -o scjoint_s1.log
#BSUB -e scjoint_s1.err


module load conda3/202210 

conda activate seurat4

R CMD BATCH 00_scjoint_s1.R

conda deactivate

