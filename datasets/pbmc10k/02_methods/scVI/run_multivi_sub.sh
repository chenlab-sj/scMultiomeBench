#!/bin/bash

#BSUB -P scvi
#BSUB -J 2_pbmc10k_control
#BSUB -q standard
#BSUB -R rusage[mem=300000]
#BSUB -n 2
#BSUB -o 2_pbmc10k_control_scvi.log
#BSUB -e 2_pbmc10k_control_scvi.err


module load conda3/202210
module load gcc/12.2.0

conda activate scvi-env
python run_multivi.py

conda deactivate
