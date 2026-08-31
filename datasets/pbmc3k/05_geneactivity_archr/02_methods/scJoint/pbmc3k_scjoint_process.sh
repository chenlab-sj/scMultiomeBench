#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P scjoint
#BSUB -J pbmc3k_scjoint_archr
#BSUB -q dgx
#BSUB -R rusage[mem=100000]
#BSUB -gpu "num=1"
#BSUB -o pbmc3k_scjointgpu_archr.log
#BSUB -e pbmc3k_scjointgpu_archr.err


module load conda3/202210

source activate scjoint
module load gcc/13.1.0
# cd into the scJoint code root so `import process_db` resolves; the process script is the ArchR variant.
cd ${TOOLS_ROOT}/scJoint/
PYTHONPATH=${TOOLS_ROOT}/scJoint:$PYTHONPATH python ${PROJECT_ROOT}/pbmc/pbmc3k/benchmark/geneactivity_archr/scJoint/pbmc3k_scjoint_process.py

conda deactivate
