#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P benchmark
#BSUB -J scjoint_s2
#BSUB -q rhel8_gpu
#BSUB -R rusage[mem=100000]
#BSUB -gpu "num=1"
#BSUB -o scjoint_s2.log
#BSUB -e scjoint_s2.err


module load conda3/202210

source activate scjoint
module load gcc/13.1.0
cp scjoint_s2.py ${TOOLS_ROOT}/scJoint/
cd ${TOOLS_ROOT}/scJoint/
python scjoint_s2.py

conda deactivate
