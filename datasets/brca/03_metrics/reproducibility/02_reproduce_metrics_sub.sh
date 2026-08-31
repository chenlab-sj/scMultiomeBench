#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J brca_fig4_reproduce
#BSUB -q large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o fig4_reproduce.log
#BSUB -e fig4_reproduce.err

# BRCA Fig4A/B reproducibility: KNN-label std (Fig4A) + pairwise-NMI (Fig4B) over 14 methods x 3 reps.
# benchmark_env (sklearn + scanpy + benchmark_fun). The 42 louvain clusterings are the slow part.
# plot_fig4.py needs seaborn/matplotlib (same env). Run plot after compute, or on your Mac for PDFs.
module load conda3/202303
conda activate benchmark_env

ROOT=${PROJECT_ROOT}
cd "${ROOT}/BRCA/benchmark/fig4"
export BENCHMARK_FUN_DIR=${BENCHMARK_FUN_DIR}
export REPRO_OUT="."
export LABEL=label.csv

python 02_reproduce_metrics.py
python plot_fig4.py

conda deactivate
echo "DONE: reproducibility.csv + fig4a_reproducibility.pdf + fig4b_nmi.pdf"
