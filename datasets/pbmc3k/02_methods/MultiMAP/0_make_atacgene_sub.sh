#!/bin/bash

#BSUB -P multimap
#BSUB -J 0_make_atacgene
#BSUB -q rhel8_compbio
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o 0_make_atacgene.log
#BSUB -e 0_make_atacgene.err


module load conda3/202210 
conda activate seurat4

Rscript 0_make_atacgene.R

conda deactivate
