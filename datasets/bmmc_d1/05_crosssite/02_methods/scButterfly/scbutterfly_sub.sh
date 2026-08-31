#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J bmmc_run_scbutterfly
#BSUB -q dgx
#BSUB -R rusage[mem=300000]
#BSUB -gpu "num=1"
# 300G (was 150G -> TERM_MEMLIMIT): scButterfly's ATAC TFIDF densifies the 39439x183117 peak matrix
# (~58G dense, ~120G peak with the divide copy). If dgx caps below this, drop to fewer concatenated
# test batches or pre-select ATAC peaks.
#BSUB -o run_scbutterfly.log
#BSUB -e run_scbutterfly.err

# scButterfly rep1 on BMMC_d1 (GPU). Single run (no batch_key; scButterfly has no dataset-batch param).
# Avoid the skip-retrain pitfall: a tiny time= means it loaded a stale model -> rm scbutterfly/model and re-run.
module load conda3/202311
source activate scbutterfly-env
export PYTHONNOUSERSITE=1
cd ${PROJECT_ROOT}/BMMC_d1/scripts/crosssite/scbutterfly
python 01_run_scbutterfly.py
conda deactivate
echo "DONE: latent.csv"
