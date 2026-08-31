#!/bin/bash

#BSUB -P scBridge
#BSUB -J scBridge_s1
#BSUB -q rhel8_standard
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o scBridge_s1.log
#BSUB -e scBridge_s1.err


module load conda3/202210 
conda activate test_env

python 00_process_scbridge_h5.py
