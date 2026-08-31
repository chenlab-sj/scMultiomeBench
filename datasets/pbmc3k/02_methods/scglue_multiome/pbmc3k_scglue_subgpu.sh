#!/bin/bash

#BSUB -P scglue
#BSUB -J pbmc3k_scglue
#BSUB -q gpu
#BSUB -R rusage[mem=200000]
#BSUB -gpu "num=1"
#BSUB -o pbmc3k_scglue_gpu.log
#BSUB -e pbmc3k_scglue_gpu.err


module load conda3/201903
source activate scglue2

python pbmc3k_withpair_scglue.py

conda deactivate
