#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J rms_midas_reps
#BSUB -q dgx
#BSUB -R rusage[mem=150000]
#BSUB -gpu "num=1"
#BSUB -o run_midas_reps.log
#BSUB -e run_midas_reps.err

# MIDAS rep2 + rep3 (seeds 0, 42) -> midas/rep2,rep3/latent.csv. Serial; split into two jobs if the
# queue time limit is tight (RMS is mid-sized -> ~between pbmc3k and BRCA in train time).
module load conda3/202311
source activate midas-env
export PYTHONNOUSERSITE=1
cd ${PROJECT_ROOT}/RMS/Mast607/script/midas
for rs in "rep2 0" "rep3 42"; do
  set -- $rs; export REP=$1 SEED=$2
  echo "=== MIDAS $REP (seed $SEED) ==="; python 01_run_midas.py
done
conda deactivate
echo "DONE: midas/rep2,rep3/latent.csv"
