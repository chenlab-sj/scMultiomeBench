#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P benchmark
#BSUB -J ga_archr_adjaccu
#BSUB -q standard
#BSUB -R rusage[mem=20000]
#BSUB -n 1
#BSUB -o ga_archr_adjaccu.log
#BSUB -e ga_archr_adjaccu.err
module load conda3/202303
conda activate benchmark_env
cd ${PROJECT_ROOT}/pbmc/pbmc3k/benchmark/geneactivity_archr/metrics
python 01_adjust_accuracy.py --label label.csv --out .
conda deactivate
echo "DONE: adj_atac_predaccu.csv"
