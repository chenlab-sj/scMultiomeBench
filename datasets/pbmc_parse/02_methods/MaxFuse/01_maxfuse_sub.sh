#!/bin/bash

#BSUB -P benchmark
#BSUB -J run_maxfuse
#BSUB -q large_mem
#BSUB -R rusage[mem=100000]
#BSUB -n 4
#BSUB -o run_maxfuse.log
#BSUB -e run_maxfuse.err

# MaxFuse is CPU-only (no GPU). If 'standard' isn't a valid queue here, swap it
# for your usual CPU queue (e.g. large_mem).
module load conda3/202311

source activate maxfuse-env

# Isolate from ~/.local user-site so the env's packages aren't shadowed.
export PYTHONNOUSERSITE=1

python 01_run_maxfuse.py

conda deactivate
