#!/bin/bash

#BSUB -P scvi
#BSUB -J run_scvi
#BSUB -q rhel8_standard
#BSUB -R rusage[mem=100000]
#BSUB -n 2
#BSUB -o run_scvi.log
#BSUB -e run_scvi.err


module load conda3/202210
module load gcc/12.2.0

conda activate scvi-env

python 01_run_scvi.py 

conda deactivate
