#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P scjoint
#BSUB -J scjoint_s0_archr
#BSUB -q compbio
#BSUB -R rusage[mem=120000]
#BSUB -n 2
#BSUB -o 0_make_scjoint_h5_archr.log
#BSUB -e 0_make_scjoint_h5_archr.err

module load conda3/202210
conda activate seurat4
cd ${PROJECT_ROOT}/pbmc/pbmc3k/benchmark/geneactivity_archr/scJoint
R CMD BATCH 00_make_scjoint_h5.R
