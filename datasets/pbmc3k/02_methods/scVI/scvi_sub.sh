#!/bin/bash

#BSUB -P scvi
#BSUB -J pbmc3k_scvi
#BSUB -q rhel8_standard
#BSUB -R rusage[mem=100000]
#BSUB -n 2
#BSUB -o pbmc3k_scvi.log
#BSUB -e pbmc3k_scvi.err


module load conda3/202210
module load gcc/12.2.0

conda activate scvi-env

python pbmc3k_scvi.py 

conda deactivate
