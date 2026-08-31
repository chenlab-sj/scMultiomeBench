#!/bin/bash

#BSUB -P scglue
#BSUB -J 2_pbmc10k_control_scglue
#BSUB -q gpu
#BSUB -R rusage[mem=200000]
#BSUB -gpu "num=1"
#BSUB -o 2_pbmc10k_control_scglue_gpu.log
#BSUB -e 2_pbmc10k_control_scglue_gpu.err


module load conda3/201903
source activate scglue2

python 01_run_scglue.py

conda deactivate
