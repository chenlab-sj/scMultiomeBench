#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J rms_scjoint_reps
#BSUB -q dgx
#BSUB -R rusage[mem=300000]
#BSUB -gpu "num=1"
#BSUB -o run_scjoint_reps.log
#BSUB -e run_scjoint_reps.err

# scJoint rep2 + rep3 (seeds 7, 42). Mirrors BRCA's 02_scjoint_s3_sub.sh: the s1/s2 prep (the npz inputs
# in ${TOOLS_ROOT}/scJoint/RMS/Mast607A_TB19_22652/) is SHARED + already done, so a rep is just
# re-running main.py with a different seed. 02_config.py now reads self.seed from env SEED. Per rep: copy
# the parameterized 02_config.py into the cloned scJoint repo, run main.py -> output/ embeddings, copy
# output/ into the rep dir. rep1 = the existing output/ (untouched). The reproduce step reads each rep's
# output/{rna_scjoint_embeddings.txt, atac_gene_scjoint_embeddings.txt} (scJoint has no single latent.csv).
module load conda3/202210
source activate scjoint
module load gcc/13.1.0

SCJOINT=${TOOLS_ROOT}/scJoint
REPBASE=${PROJECT_ROOT}/RMS/Mast607/script/scJoint

for rs in "rep2 7" "rep3 42"; do
  set -- $rs; export REP=$1 SEED=$2
  echo "=== scJoint $REP (seed $SEED) ==="
  outdir="$REPBASE/$REP"
  mkdir -p "$outdir"
  cp "$REPBASE/config.py" "$SCJOINT/config.py"   # parameterized: reads SEED from env at import
  cd "$SCJOINT"
  python main.py
  cp -r output/ "$outdir/"
done

conda deactivate
echo "DONE: scJoint rep2,rep3 output/ embeddings"
