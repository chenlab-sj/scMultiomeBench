#!/bin/bash

#BSUB -P MinNet
#BSUB -J 1_MinNet_train_MinNet
#BSUB -q gpu
#BSUB -R rusage[mem=200000]
#BSUB -gpu "num=1"
#BSUB -o 1_MinNet_train_MinNet.log
#BSUB -e 1_MinNet_train_MinNet.err

module load conda3/202105

source activate MinNet

python 1_traindata_process.py