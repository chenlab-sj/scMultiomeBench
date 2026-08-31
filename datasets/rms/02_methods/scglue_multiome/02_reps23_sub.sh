#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P scglue
#BSUB -J rms_scglue_withpair_reps
#BSUB -q dgx
#BSUB -R rusage[mem=150000]
#BSUB -gpu "num=1"
#BSUB -o run_scglue_withpair_reps.log
#BSUB -e run_scglue_withpair_reps.err

# scGLUE (with pair / PairedSCGLUEModel) rep2 + rep3 (seeds 7, 42) -> out_dir/rep2,rep3/scglue_latent.csv.
# rep1 = the existing out_dir/scglue_latent.csv (untouched). 7/42 avoid the original seed 420.
# GPU -> dgx (matches 01_run_scglue_sub.sh). REP routes the latent AND all intermediates
# (rna/atac-pp.h5ad, guidance.graphml.gz, glue.dill, glue/ training dir) under out_dir/<REP>
# so the two reps never clobber each other or rep1.
module load conda3/201903
source activate scglue2
cd ${PROJECT_ROOT}/RMS/Mast607/script/scglue_withpair
for rs in "rep2 7" "rep3 42"; do
  set -- $rs; export REP=$1 SEED=$2
  echo "=== scglue_withpair $REP (seed $SEED) ==="; python 01_run_scglue_multiome.py
done
conda deactivate
echo "DONE: out_dir/rep2,rep3/scglue_latent.csv"
