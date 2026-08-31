#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P bindsc
#BSUB -J pbmc3k_bindsc_archr
#BSUB -q large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o pbmc3k_bindsc_archr.log
#BSUB -e pbmc3k_bindsc_archr.err


module load conda3/202210

conda activate bindsc

cd ${PROJECT_ROOT}/pbmc/pbmc3k/benchmark/geneactivity_archr/bindsc

R CMD BATCH pbmc3k_bindsc.R

conda deactivate
