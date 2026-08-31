#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J scBridge_s1
#BSUB -q gpu_short
#BSUB -R rusage[mem=300000]
#BSUB -gpu "num=1"
#BSUB -o scBridge_s1.log
#BSUB -e scBridge_s1.err


module load conda3/202210

conda activate scBridge
module load gcc/13.1.0-rhel7

cd ${TOOLS_ROOT}/scBridge/
python main.py --data_path="${PROJECT_ROOT}/BMMC_d1/s1d1_paired/scBridge/" --source_data="scBridge_rna.h5ad" --target_data="scBridge_atac.h5ad" --umap_plot

conda deactivate
