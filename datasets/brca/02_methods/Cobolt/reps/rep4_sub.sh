#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P cobolt
#BSUB -J brca_cobolt_rep4
#BSUB -q gpu
#BSUB -R rusage[mem=100000]
#BSUB -gpu "num=1"
#BSUB -o run_cobolt_rep4.log
#BSUB -e run_cobolt_rep4.err

# BRCA Cobolt rep4 (seed 100), ~5h. Writes repo/cobolt/rep4/{latent,cobolt_latent}.csv. Run in parallel
# with rep5_sub.sh (each its own GPU) so neither hits a serial-walltime wall.
module load conda3/202210
conda activate cobolt
cd ${PROJECT_ROOT}/BRCA/HT243B1-S1H4/cobolt
mkdir -p rep4
export REP=rep4 SEED=100
python run_cobolt_reps.py
conda deactivate
