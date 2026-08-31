#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J rms_scbutterfly_reps
#BSUB -q dgx
#BSUB -R rusage[mem=150000]
#BSUB -gpu "num=1"
#BSUB -o run_scbutterfly_reps.log
#BSUB -e run_scbutterfly_reps.err

# scButterfly rep2 + rep3 (seeds 0, 42). gpu_interactive (NOT gpu_short -> avoids the TERM_RUNLIMIT that
# killed pbmc3k rep3). Each rep trains into its OWN repN/model/. Split into two jobs if time is tight.
module load conda3/202311
source activate scbutterfly-env
export PYTHONNOUSERSITE=1
cd ${PROJECT_ROOT}/RMS/Mast607/script/scbutterfly
for rs in "rep2 0" "rep3 42"; do
  set -- $rs; export REP=$1 SEED=$2
  echo "=== scButterfly $REP (seed $SEED) ==="; python 01_run_scbutterfly.py
done
conda deactivate
echo "DONE: scbutterfly/rep2,rep3/latent.csv"
