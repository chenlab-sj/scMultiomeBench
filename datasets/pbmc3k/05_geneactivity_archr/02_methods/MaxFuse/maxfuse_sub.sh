#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J run_maxfuse_archr
#BSUB -q large_mem
#BSUB -R rusage[mem=100000]
#BSUB -n 4
#BSUB -o run_maxfuse_archr.log
#BSUB -e run_maxfuse_archr.err

# MaxFuse is CPU-only (no GPU). If 'standard' isn't a valid queue here, swap it
# for your usual CPU queue (e.g. large_mem).
cd ${PROJECT_ROOT}/pbmc/pbmc3k/benchmark/geneactivity_archr/maxfuse

module load conda3/202311

source activate maxfuse-env

# Isolate from ~/.local user-site so the env's packages aren't shadowed.
export PYTHONNOUSERSITE=1

python 01_run_maxfuse.py

conda deactivate
