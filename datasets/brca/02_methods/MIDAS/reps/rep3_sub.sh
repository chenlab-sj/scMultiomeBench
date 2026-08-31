#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J brca_midas_rep3
#BSUB -q dgx
#BSUB -R rusage[mem=200000]
#BSUB -gpu "num=1"
#BSUB -o run_midas_rep3.log
#BSUB -e run_midas_rep3.err

# MIDAS rep3 (seed 42) -> midas/rep3/latent.csv. Single-rep parallel job, run alongside the
# serial 02_reps23_sub.sh (which is busy training rep2).
# IMPORTANT: the serial 02_reps23_sub.sh runs rep2 THEN rep3. Once midas/rep2/latent.csv exists, bkill
# the serial reps job so it does NOT also start rep3 and collide with this job over midas/rep3/.
module load conda3/202311
source activate midas-env
export PYTHONNOUSERSITE=1

cd ${PROJECT_ROOT}/BRCA/HT243B1-S1H4/midas

export REP=rep3 SEED=42
python 01_run_midas.py

conda deactivate
echo "DONE: midas/rep3/latent.csv"
