#!/bin/bash

#BSUB -P benchmark
#BSUB -J run_mira_unp
#BSUB -q gpu_interactive
#BSUB -R rusage[mem=100000]
#BSUB -gpu "num=1"
#BSUB -o run_mira_unpaired.log
#BSUB -e run_mira_unpaired.err

module load conda3/202311

source activate mira-env

# Isolate from ~/.local user-site so the env's mira-multiome isn't shadowed.
export PYTHONNOUSERSITE=1

python run_mira_unpaired.py

conda deactivate
