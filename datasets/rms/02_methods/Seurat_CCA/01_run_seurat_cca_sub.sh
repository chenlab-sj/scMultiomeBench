#!/bin/bash

#BSUB -P Seurat3
#BSUB -J RMS_Mast607A_seurat3
#BSUB -q rhel8_compbio
#BSUB -R rusage[mem=100000]
#BSUB -n 2
#BSUB -o RMS_Mast607A_seurat3.log
#BSUB -e RMS_Mast607A_seurat3.err


module load conda3/202210 

conda activate seurat4

R CMD BATCH 01_run_seurat_cca.R

conda deactivate
