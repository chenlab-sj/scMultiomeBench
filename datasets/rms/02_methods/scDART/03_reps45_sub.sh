#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J rms_scdart_reps2
#BSUB -q dgx
#BSUB -gpu "num=1"
#BSUB -R rusage[mem=200000]
#BSUB -o run_scdart_reps2.log
#BSUB -e run_scdart_reps2.err

# scDART EXTRA seeds (rep4=13, rep5=99) to inspect scDART's spread across more runs.
# Seeds 13/99 avoid the ones already used: 420 (base/rep1), 7 (rep2), 42 (rep3), and the
# originals' hardcoded 0/1/1234. Same cloned-repo mechanism as 02_reps23_sub.sh: copy the
# parameterized 01_run_scdart.py into the scDART Examples/ dir and run it there; the script writes
# to .../scDART/<REP>/latent.csv. rep1/rep2/rep3 are untouched.
module load conda3/202402
conda activate scdart

REPO_SCRIPT=${PROJECT_ROOT}/RMS/Mast607/script/scDART/run_scdart.py
EXAMPLES=${TOOLS_ROOT}/scDART/scDART/Examples

for rs in "rep4 13" "rep5 99"; do
  set -- $rs; export REP=$1 SEED=$2
  echo "=== scDART $REP (seed $SEED) ==="
  cp "$REPO_SCRIPT" "$EXAMPLES/run_scdart.py"
  cd "$EXAMPLES"
  python 01_run_scdart.py
done

conda deactivate
echo "DONE: scDART/rep4,rep5/latent.csv (seeds 13, 99)"
