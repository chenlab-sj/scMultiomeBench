#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J rms_portal_reps45
#BSUB -q dgx
#BSUB -gpu "num=1"
#BSUB -R rusage[mem=150000]
#BSUB -o run_portal_reps45.log
#BSUB -e run_portal_reps45.err

# Portal rep4 + rep5 (seeds 100, 200). Portal is run from the cloned Portal repo dir, matching
# 02_reps23_sub.sh. The parameterized 01_run_portal.py writes to /home/.../portal/<REP>/lat_df.csv; this
# script then copies each latent back into this benchmark repo for 01_prep_latents.py staging.
module load conda3/202402
conda activate portal

REPO_DIR=${PROJECT_ROOT}/RMS/Mast607/script/portal
REPO_SCRIPT="${REPO_DIR}/run_portal.py"
PORTAL_REPO=${TOOLS_ROOT}/portal/Portal
PORTAL_OUT=${DATA_ROOT}/RMS/Mast607A_TB19_22652/portal

for rs in "rep4 100" "rep5 200"; do
  set -- $rs
  export REP=$1 SEED=$2
  echo "=== portal $REP (seed $SEED) ==="
  cp "$REPO_SCRIPT" "$PORTAL_REPO/run_portal.py"
  cd "$PORTAL_REPO" || exit 1
  python 01_run_portal.py

  mkdir -p "${REPO_DIR}/${REP}"
  cp "${PORTAL_OUT}/${REP}/lat_df.csv" "${REPO_DIR}/${REP}/lat_df.csv"
done

conda deactivate
echo "DONE: RMS/Mast607/script/portal/rep4,rep5/lat_df.csv"
