#!/bin/bash

#BSUB -P rliger
#BSUB -J pbmc10k_rliger_control
#BSUB -q standard
#BSUB -R rusage[mem=200000]
#BSUB -n 1
#BSUB -o 2_pbmc10k_control_liger.log
#BSUB -e 2_pbmc10k_control_liger.err

module load conda3/202210

conda activate liger
R CMD BATCH run_liger.R

conda deactivate