#!/bin/bash

#BSUB -P benchmark
#BSUB -J run_midas
#BSUB -q dgx
#BSUB -R rusage[mem=100000]
#BSUB -gpu "num=1"
#BSUB -o run_midas.log
#BSUB -e run_midas.err

module load conda3/202311

source activate midas-env

export PYTHONNOUSERSITE=1

python 01_run_midas.py

conda deactivate
