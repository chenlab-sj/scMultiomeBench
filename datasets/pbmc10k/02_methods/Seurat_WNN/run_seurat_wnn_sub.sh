#!/bin/bash

#BSUB -P seurat4
#BSUB -J 1_3_pbmc10k_control_seurat4
#BSUB -q standard
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o 1_3_pbmc10k_control_seurat4.log
#BSUB -e 1_3_pbmc10k_control_seurat4.err


module load conda3/202210 
conda activate seurat4

R CMD BATCH run_seurat_wnn.R

conda deactivate
