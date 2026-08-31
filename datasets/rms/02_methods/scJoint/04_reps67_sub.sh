#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J rms_scjoint_reps67
#BSUB -q dgx
#BSUB -R rusage[mem=300000]
#BSUB -gpu "num=1"
#BSUB -o run_scjoint_reps67.log
#BSUB -e run_scjoint_reps67.err

# scJoint rep6 + rep7 (seeds 13, 99). The shared scJoint npz inputs are already
# prepared in ${TOOLS_ROOT}/scJoint/RMS/Mast607A_TB19_22652/.
set -euo pipefail

module load conda3/202210
source activate scjoint
module load gcc/13.1.0

SCJOINT=${TOOLS_ROOT}/scJoint
REPBASE=${PROJECT_ROOT}/RMS/Mast607/script/scJoint

for rs in "rep6 13" "rep7 99"; do
  set -- $rs; export REP=$1 SEED=$2
  echo "=== scJoint $REP (seed $SEED) ==="
  outdir="$REPBASE/$REP/output"
  mkdir -p "$outdir"
  cp "$REPBASE/config.py" "$SCJOINT/config.py"
  cd "$SCJOINT"
  python main.py
  cp output/* "$outdir/"
done

conda deactivate
echo "DONE: scJoint rep6,rep7 output/ embeddings"
