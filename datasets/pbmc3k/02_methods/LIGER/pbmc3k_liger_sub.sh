#!/bin/bash

#BSUB -P rliger
#BSUB -J pbmc3k_liger
#BSUB -q standard
#BSUB -R rusage[mem=200000]
#BSUB -n 1
#BSUB -o pbmc3k_liger.log
#BSUB -e pbmc3k_liger.err

module load conda3/202210

conda activate liger
R CMD BATCH pbmc3k_liger.R

conda deactivate