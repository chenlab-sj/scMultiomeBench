#!/bin/bash

#BSUB -P multimap
#BSUB -J 2_pbmc10k_control_multimap
#BSUB -q standard
#BSUB -R rusage[mem=100000]
#BSUB -n 1
#BSUB -o 2_pbmc10k_control_multimap.log
#BSUB -e 2_pbmc10k_control_multimap.err


module load conda3/202210 

source activate multimap
python 01_run_multimap.py

conda deactivate

