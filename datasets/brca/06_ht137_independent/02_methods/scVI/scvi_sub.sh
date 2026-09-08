#!/bin/bash

#BSUB -P benchmark
#BSUB -J run_scVI
#BSUB -q gpu
#BSUB -R rusage[mem=100000]
#BSUB -gpu "num=1"
#BSUB -o run_scVIgpu.log
#BSUB -e run_scVIgpu.err

module unload conda3/202210
module load conda3/202311

conda activate scvi-env

python run_scvi.py 

conda deactivate
