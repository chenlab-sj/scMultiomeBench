#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J rms_portal_reps
#BSUB -q dgx
#BSUB -gpu "num=1"
#BSUB -R rusage[mem=150000]
#BSUB -o run_portal_reps.log
#BSUB -e run_portal_reps.err

# Portal rep2 + rep3 (seeds 7, 42). Portal is run from the CLONED Portal repo dir: the original
# run_portal_sub.sh does `cp 01_run_portal.py ${TOOLS_ROOT}/portal/Portal` then `cd` there and runs.
# So, exactly like the original sub, copy the (parameterized) 01_run_portal.py into the repo and run
# it there. rep1 = the existing portal/lat_df.csv (untouched). The parameterized script writes to
# .../portal/<REP>/lat_df.csv. 7/42 avoid the original seed 420. GPU (-gpu) -> dgx.
module load conda3/202402
conda activate portal

REPO_SCRIPT=${PROJECT_ROOT}/RMS/Mast607/script/portal/run_portal.py
PORTAL_REPO=${TOOLS_ROOT}/portal/Portal

for rs in "rep2 7" "rep3 42"; do
  set -- $rs; export REP=$1 SEED=$2
  echo "=== portal $REP (seed $SEED) ==="
  cp "$REPO_SCRIPT" "$PORTAL_REPO/run_portal.py"
  cd "$PORTAL_REPO"
  python 01_run_portal.py
done

conda deactivate
echo "DONE: ${DATA_ROOT}/RMS/Mast607A_TB19_22652/portal/rep2,rep3/lat_df.csv"
