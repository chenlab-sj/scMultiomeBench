#!/bin/bash

#BSUB -P cobolt
#BSUB -J run_cobolt
#BSUB -q dgx
#BSUB -R rusage[mem=100000]
#BSUB -gpu "num=1"
#BSUB -o run_coboltgpu.log
#BSUB -e run_coboltgpu.err


module load conda3/202210

conda activate cobolt
python 01_run_cobolt.py

conda deactivate
