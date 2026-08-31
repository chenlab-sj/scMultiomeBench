#!/bin/bash

#BSUB -P benchmark
#BSUB -J bindsc
#BSUB -q standard
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o bindsc.log
#BSUB -e bindsc.err


module load conda3/202210 

conda activate bindsc

R CMD BATCH bindsc.R

conda deactivate

