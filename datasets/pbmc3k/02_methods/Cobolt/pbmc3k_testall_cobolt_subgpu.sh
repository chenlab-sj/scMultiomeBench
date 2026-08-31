#!/bin/bash

#BSUB -P cobolt
#BSUB -J pbmc3k_test_all_cobolt
#BSUB -q rhel8_gpu
#BSUB -R rusage[mem=100000]
#BSUB -gpu "num=1"
#BSUB -o pbmc3k_test_all_coboltgpu.log
#BSUB -e pbmc3k_test_all_coboltgpu.err


module load conda3/202210

conda activate cobolt
python pbmc3k_test_all_cobolt.py

conda deactivate
