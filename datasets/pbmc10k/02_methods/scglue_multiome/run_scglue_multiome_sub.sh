#!/bin/bash

#BSUB -P scglue
#BSUB -J 4_pbmc10k_withpair_control_scglue
#BSUB -q gpu
#BSUB -R rusage[mem=200000]
#BSUB -gpu "num=1"
#BSUB -o 4_pbmc10k_withpair_control_scglue_gpu.log
#BSUB -e 4_pbmc10k_withpair_control_scglue_gpu.err


module load conda3/201903
source activate scglue2

python run_scglue_multiome.py

conda deactivate
