#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J scjoint_s3
#BSUB -q gpu_short
#BSUB -R rusage[mem=300000]
#BSUB -gpu "num=1"
#BSUB -o scjointgpu_s3.log
#BSUB -e scjointgpu_s3.err


module load conda3/202210

source activate scjoint
module load gcc/13.1.0

outdir=${LS_SUBCWD:-$PWD}
cd "$outdir" || exit 1
cp 02_config.py ${TOOLS_ROOT}/scJoint/ || exit 1
cd ${TOOLS_ROOT}/scJoint || exit 1

# ${TOOLS_ROOT}/scJoint is SHARED by every scJoint run in this repo (it still holds pbmc3k_*
# leftovers). Delete the two embedding files we are about to regenerate, so that if main.py fails
# s4 errors out instead of silently assembling a latent from the previous run's dataset -- which
# is exactly how the 8291-vs-8442 assert happened.
# NOTE: for the same reason, crosssite and s1d1_paired s3 jobs must NOT run at the same time.
rm -f output/rna_scjoint_embeddings.txt output/atac_gene_scjoint_embeddings.txt

python main.py
rc=$?
if [ $rc -ne 0 ]; then
  echo "ERROR: scJoint main.py failed (rc=$rc)"
  conda deactivate
  exit $rc
fi

rm -rf "$outdir/output"
cp -r output/ "$outdir" || exit 1

python "$outdir/scjoint_s4_assemble.py"
rc=$?

conda deactivate
# Propagate the real status: previously this script ended on `conda deactivate`, so LSF reported
# "Successfully completed" even when main.py and s4 had both crashed.
[ $rc -eq 0 ] && echo "DONE: scJoint_latent.csv"
exit $rc

