#!/bin/bash

#BSUB -P MinNet
#BSUB -J 1_2_MinNet_test
#BSUB -q gpu
#BSUB -R rusage[mem=200000]
#BSUB -gpu "num=1"
#BSUB -o 1_2_MinNet_test.log
#BSUB -e 1_2_MinNet_test.err

module load conda3/202105

source activate MinNet

python 1_2_test_dataprocess.py