#!/bin/bash

#BSUB -P bindsc
#BSUB -J run_bindsc
#BSUB -q rhel8_large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o run_bindsc.log
#BSUB -e run_bindsc.err


module load conda3/202210 

conda activate bindsc

R CMD BATCH 01_run_bindsc.R

conda deactivate
