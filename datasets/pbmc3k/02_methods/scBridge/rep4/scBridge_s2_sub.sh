#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P scBridge
#BSUB -J pbmc3k_scb_s2_rep4
#BSUB -q large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o scBridge_s2.log
#BSUB -e scBridge_s2.err

# scBridge pbmc3k reproducibility rep4, step 2: read back the embedding (tmp pkl from step 1) and write
# pbmc3k/latent.csv in the LSA scBridge dir, then copy it into THIS rep dir (the repo archive that
# reproduce_extrareps.py reads). Submit from this dir AFTER 01_scBridge_s1_sub.sh completes:
#   cd scripts/scBridge/rep4 && bsub < 02_scBridge_s2_sub.sh
module load conda3/202210
conda activate scBridge
module load gcc/13.1.0-rhel7

outdir=$(pwd)
cd ${TOOLS_ROOT}/scBridge
python 2_save_scBridge_result.py                 # -> pbmc3k/latent.csv, pbmc3k/label_pred.csv
cp pbmc3k/latent.csv     "$outdir/latent.csv"
cp pbmc3k/label_pred.csv "$outdir/label_pred.csv"
conda deactivate
