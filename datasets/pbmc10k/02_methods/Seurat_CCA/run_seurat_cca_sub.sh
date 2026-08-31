#!/bin/bash

#BSUB -P Seuratv3
#BSUB -J 2_pbmc10k_control_seurat3
#BSUB -q standard
#BSUB -R rusage[mem=120000]
#BSUB -n 2
#BSUB -o 2_pbmc10k_control_seurat3.log
#BSUB -e 2_pbmc10k_control_seurat3.err


module load conda3/202210 

conda activate seurat4

R CMD BATCH run_seurat_cca.R

conda deactivate
