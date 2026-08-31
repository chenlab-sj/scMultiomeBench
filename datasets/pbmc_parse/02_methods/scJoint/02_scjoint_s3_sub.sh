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

outdir=$(pwd)
cp 02_config.py ${TOOLS_ROOT}/scJoint/
cd ${TOOLS_ROOT}/scJoint

python main.py

cp -r output/ $outdir
conda deactivate

