#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J bmmc_run_maxfuse
#BSUB -q large_mem
#BSUB -R rusage[mem=150000]
#BSUB -n 4
#BSUB -o run_maxfuse.log
#BSUB -e run_maxfuse.err

# MaxFuse Fusor on BMMC_d1 (CPU). Reuses the BMMC gene-activity / RNA / LSI h5ads built by
# 00_file_prep.R (test_rna.h5ad, test_atac_gene.h5ad, atac_embed.h5ad) -> no separate R prep step.
module load conda3/202311
source activate maxfuse-env
export PYTHONNOUSERSITE=1
cd ${PROJECT_ROOT}/BMMC_d1/scripts/maxfuse
python 01_run_maxfuse.py
conda deactivate
echo "DONE: latent.csv"
