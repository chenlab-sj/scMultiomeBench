#!/bin/bash

#BSUB -P benchmark
#BSUB -J run_scbutterfly
#BSUB -q gpu_short
#BSUB -R rusage[mem=100000]
#BSUB -gpu "num=1"
#BSUB -o run_scbutterfly.log
#BSUB -e run_scbutterfly.err

module load conda3/202311

source activate scbutterfly-env

export PYTHONNOUSERSITE=1

python 01_run_scbutterfly.py

conda deactivate
