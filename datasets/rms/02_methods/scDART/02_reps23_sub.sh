#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J rms_scdart_reps
#BSUB -q dgx
#BSUB -gpu "num=1"
#BSUB -R rusage[mem=200000]
#BSUB -o run_scdart_reps.log
#BSUB -e run_scdart_reps.err

# scDART rep2 + rep3 (seeds 7, 42). scDART is a CLONED repo run from its Examples/ dir (01_run_scdart.py
# does sys.path.append('../') so the scDART package only imports from there). So, exactly like the
# original run_scdart_sub.sh, copy the (parameterized) 01_run_scdart.py into Examples/ and run it there.
# rep1 = the existing latent.csv (untouched). The parameterized script writes to .../scDART/<REP>/.
module load conda3/202402
conda activate scdart

REPO_SCRIPT=${PROJECT_ROOT}/RMS/Mast607/script/scDART/run_scdart.py
EXAMPLES=${TOOLS_ROOT}/scDART/scDART/Examples

for rs in "rep2 7" "rep3 42"; do
  set -- $rs; export REP=$1 SEED=$2
  echo "=== scDART $REP (seed $SEED) ==="
  cp "$REPO_SCRIPT" "$EXAMPLES/run_scdart.py"
  cd "$EXAMPLES"
  python 01_run_scdart.py
done

conda deactivate
echo "DONE: ${DATA_ROOT}/RMS/Mast607A_TB19_22652/scDART/rep2,rep3/latent.csv"
