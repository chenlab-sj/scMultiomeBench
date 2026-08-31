#!/bin/bash

#BSUB -P scjoint
#BSUB -J pbmc10k_control_scjoint
#BSUB -q rhel8_gpu
#BSUB -R rusage[mem=100000]
#BSUB -gpu "num=1"
#BSUB -o 2_pbmc10k_control_scjointgpu.log
#BSUB -e 2_pbmc10k_control_scjointgpu.err


module load conda3/202210

source activate scjoint
module load gcc/13.1.0
python 01_scjoint_data_parsing.py

conda deactivate
