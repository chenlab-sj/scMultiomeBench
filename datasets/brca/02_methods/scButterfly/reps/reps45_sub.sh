#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J brca_scbutterfly_reps45
#BSUB -q gpu_interactive
#BSUB -R rusage[mem=200000]
#BSUB -gpu "num=1"
#BSUB -o run_scbutterfly_reps45.log
#BSUB -e run_scbutterfly_reps45.err

# scButterfly rep4 + rep5 (seeds 100, 200), FORCE_SEED=1 overrides scButterfly internal seed 19193 -> scbutterfly/rep4/latent.csv, scbutterfly/rep5/latent.csv.
# Two full VAE trainings run SERIALLY. gpu_interactive for the longer wall-time; split into two
# jobs if you want them parallel. Each rep retrains (its own scbutterfly/repN/model/).
module load conda3/202311
source activate scbutterfly-env
export PYTHONNOUSERSITE=1
export FORCE_SEED=1

cd ${PROJECT_ROOT}/BRCA/HT243B1-S1H4/scbutterfly

for rs in "rep4 100" "rep5 200"; do
  set -- $rs; export REP=$1 SEED=$2
  echo "=== scButterfly $REP (seed $SEED) ==="
  python 01_run_scbutterfly.py
done

conda deactivate
echo "DONE: scbutterfly/rep2/latent.csv + scbutterfly/rep3/latent.csv"
