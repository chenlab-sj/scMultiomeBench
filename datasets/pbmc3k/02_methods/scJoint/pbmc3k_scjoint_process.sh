#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P scjoint
#BSUB -J pbmc3k_scjoint
#BSUB -q rhel8_gpu
#BSUB -R rusage[mem=100000]
#BSUB -gpu "num=1"
#BSUB -o pbmc3k_scjointgpu.log
#BSUB -e pbmc3k_scjointgpu.err


module load conda3/202210

source activate scjoint
module load gcc/13.1.0
cd ${TOOLS_ROOT}/scJoint/
python pbmc3k/pbmc3k_scjoint_process.py

conda deactivate
