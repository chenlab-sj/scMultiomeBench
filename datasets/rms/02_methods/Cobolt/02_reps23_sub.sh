#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J rms_cobolt_reps
#BSUB -q dgx
#BSUB -gpu "num=1"
#BSUB -R rusage[mem=150000]
#BSUB -o run_cobolt_reps.log
#BSUB -e run_cobolt_reps.err

# cobolt rep2 + rep3 (seeds 7, 42) -> res_Mast607A/rep2,rep3/latent.csv. rep1 = existing
# res_Mast607A/latent.csv (untouched). GPU -> dgx. 7/42 avoid the original seed 420.
module load conda3/202210
conda activate cobolt
cd ${PROJECT_ROOT}/RMS/Mast607/script/cobolt
for rs in "rep2 7" "rep3 42"; do
  set -- $rs; export REP=$1 SEED=$2
  echo "=== cobolt $REP (seed $SEED) ==="; python 01_run_cobolt.py
done
conda deactivate
echo "DONE: res_Mast607A/rep2,rep3/latent.csv"
