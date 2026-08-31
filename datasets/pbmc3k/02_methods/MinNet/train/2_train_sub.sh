#!/bin/bash

#BSUB -P MinNet
#BSUB -J 2_MinNet_train_MinNet
#BSUB -q gpu
#BSUB -R rusage[mem=200000]
#BSUB -gpu "num=1"
#BSUB -o 2_MinNet_train_MinNet.log
#BSUB -e 2_MinNet_train_MinNet.err

module load conda3/202105

source activate MinNet

python 2_training_main.py