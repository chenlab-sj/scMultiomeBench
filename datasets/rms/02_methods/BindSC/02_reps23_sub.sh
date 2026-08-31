#!/bin/bash

#BSUB -P bindsc
#BSUB -J run_bindsc_reps
#BSUB -q large_mem
#BSUB -R rusage[mem=150000]
#BSUB -n 1
#BSUB -o run_bindsc_reps.log
#BSUB -e run_bindsc_reps.err


module load conda3/202210

conda activate bindsc

for rs in "rep2 7" "rep3 42"; do
  set -- $rs
  export REP=$1 SEED=$2
  echo "=== bindsc $REP (seed $SEED) ==="
  R CMD BATCH 01_run_bindsc.R run_bindsc_${REP}.Rout
done

conda deactivate
