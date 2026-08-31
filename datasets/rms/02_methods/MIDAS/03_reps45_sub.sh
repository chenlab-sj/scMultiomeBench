#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J rms_midas_reps2
#BSUB -q dgx
#BSUB -R rusage[mem=150000]
#BSUB -gpu "num=1"
#BSUB -o run_midas_reps2.log
#BSUB -e run_midas_reps2.err

# MIDAS EXTRA seeds (rep4=13, rep5=99) -- the seed-fragility spot-check (same seeds as scDART rep4/rep5),
# to see whether MIDAS also has a collapse seed that 3 reps missed. -> midas/rep4,rep5/latent.csv.
# rep1/rep2/rep3 untouched. Serial; split into two jobs if the queue time limit is tight.
module load conda3/202311
source activate midas-env
export PYTHONNOUSERSITE=1
cd ${PROJECT_ROOT}/RMS/Mast607/script/midas
for rs in "rep4 13" "rep5 99"; do
  set -- $rs; export REP=$1 SEED=$2
  echo "=== MIDAS $REP (seed $SEED) ==="; python 01_run_midas.py
done
conda deactivate
echo "DONE: midas/rep4,rep5/latent.csv"
