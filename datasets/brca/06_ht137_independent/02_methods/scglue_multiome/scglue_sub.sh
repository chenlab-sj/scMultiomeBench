#!/bin/bash

#BSUB -P scglue
#BSUB -J scglue
#BSUB -q gpu
#BSUB -R rusage[mem=800000]
#BSUB -gpu "num=1"
#BSUB -o scglue_gpu.log
#BSUB -e scglue_gpu.err


module load conda3/201903
source activate scglue2

python withpair_scglue.py

conda deactivate

