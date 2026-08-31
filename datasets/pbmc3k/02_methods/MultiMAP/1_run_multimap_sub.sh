#!/bin/bash

#BSUB -P multimap
#BSUB -J 1_run_multimap
#BSUB -q rhel8_compbio 
#BSUB -R rusage[mem=100000]
#BSUB -n 1
#BSUB -o 1_run_multimap.log
#BSUB -e 1_run_multimap.err


module load conda3/202210 

source activate multimap
python 1_run_multimap.py

conda deactivate

