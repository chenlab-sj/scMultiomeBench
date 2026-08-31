#!/bin/bash

#BSUB -P bindsc
#BSUB -J pbmc3k_bindsc
#BSUB -q rhel8_large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o pbmc3k_bindsc.log
#BSUB -e pbmc3k_bindsc.err


module load conda3/202210 

conda activate bindsc

R CMD BATCH pbmc3k_bindsc.R

conda deactivate
