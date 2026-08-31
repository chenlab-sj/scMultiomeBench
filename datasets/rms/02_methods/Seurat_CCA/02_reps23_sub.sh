#!/bin/bash

#BSUB -P Seurat3
#BSUB -J RMS_Mast607A_seurat3_reps
#BSUB -q large_mem
#BSUB -R rusage[mem=150000]
#BSUB -n 2
#BSUB -o run_seurat_reps.log
#BSUB -e run_seurat_reps.err


module load conda3/202210

conda activate seurat4

for rs in "rep2 7" "rep3 42"; do
  set -- $rs
  export REP=$1 SEED=$2
  echo "=== seurat $REP (seed $SEED) ==="
  R CMD BATCH 01_run_seurat_cca.R RMS_Mast607A_seurat3_${REP}.Rout
done

conda deactivate
