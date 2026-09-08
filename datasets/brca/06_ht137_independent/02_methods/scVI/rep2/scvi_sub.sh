#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P benchmark
#BSUB -J scVI_rep2
#BSUB -q gpu
#BSUB -R rusage[mem=100000]
#BSUB -gpu "num=1"
#BSUB -o run_scVIgpu_rep2.log
#BSUB -e run_scVIgpu_rep2.err

module unload conda3/202210
module load conda3/202311

conda activate scvi-env
cd ${PROJECT_ROOT}/BRCA/HT137B1-S1H7/scVI/rep2

python run_scvi.py 

conda deactivate
