#!/bin/bash

#BSUB -P scglue
#BSUB -J scglue
#BSUB -q rhel8_gpu
#BSUB -R rusage[mem=800000]
#BSUB -gpu "num=1"
#BSUB -o scglue_gpu.log
#BSUB -e scglue_gpu.err


module load conda3/201903
source activate scglue2

python 01_run_scglue_multiome.py

conda deactivate

