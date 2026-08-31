#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J rms_scvi_reps
#BSUB -q large_mem
#BSUB -R rusage[mem=150000]
#BSUB -n 4
#BSUB -o run_scvi_reps.log
#BSUB -e run_scvi_reps.err

# scVI rep2 + rep3 (seeds 7, 42) -> res_Mast607A/rep2,rep3/latent.csv. rep1 = the existing
# res_Mast607A/latent.csv (untouched). 7/42 avoid the original seed 420. CPU (no GPU) -> large_mem.
module load conda3/202210
module load gcc/12.2.0
conda activate scvi-env
cd ${PROJECT_ROOT}/RMS/Mast607/script/scVI
for rs in "rep2 7" "rep3 42"; do
  set -- $rs; export REP=$1 SEED=$2
  echo "=== scVI $REP (seed $SEED) ==="; python 01_run_scvi.py
done
conda deactivate
echo "DONE: res_Mast607A/rep2,rep3/latent.csv"
