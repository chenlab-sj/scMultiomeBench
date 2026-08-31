#!/bin/bash

#BSUB -P concos
#BSUB -J pbmc10k_control_concos
#BSUB -q standard
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o pbmc10k_control_concos.log
#BSUB -e pbmc10k_control_concos.err


module load conda3/202210 
conda activate concos

R CMD BATCH run_conos.R

conda deactivate
