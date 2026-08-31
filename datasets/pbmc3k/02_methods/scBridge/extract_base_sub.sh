#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P scBridge
#BSUB -J scBridge_base_lat
#BSUB -q standard
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o extract_base.log
#BSUB -e extract_base.err

# Recover scBridge base(rep1) latent -> scBridge/latent.csv from the already-integrated h5ads.
# Fast (~2 min): just reads .obsm["Embedding"] and concats, no re-run. Same env as 2_save_scBridge.
module load conda3/202210
conda activate scBridge
module load gcc/13.1.0-rhel7

cd ${PROJECT_ROOT}/pbmc/pbmc3k/scripts/scBridge
python extract_base_latent.py

conda deactivate
