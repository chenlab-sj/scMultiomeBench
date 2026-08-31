#!/bin/bash

#BSUB -P MinNet
#BSUB -J 2_MinNet_test
#BSUB -q gpu
#BSUB -R rusage[mem=200000]
#BSUB -gpu "num=1"
#BSUB -o 2_MinNet_test.log
#BSUB -e 2_MinNet_test.err

module load conda3/202105

source activate MinNet

python 02_run_minnet_siann.py 