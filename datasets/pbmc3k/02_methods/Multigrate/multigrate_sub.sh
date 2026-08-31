#!/bin/bash

#BSUB -P benchmark
#BSUB -J run_multigrate
#BSUB -q gpu_short
#BSUB -R rusage[mem=100000]
#BSUB -gpu "num=1"
#BSUB -o run_multigrate.log
#BSUB -e run_multigrate.err

module load conda3/202311

source activate multigrate-env

export PYTHONNOUSERSITE=1

python run_multigrate.py

conda deactivate
