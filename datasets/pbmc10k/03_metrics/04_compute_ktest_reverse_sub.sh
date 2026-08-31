#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P benchmark
#BSUB -J figS1b
#BSUB -q large_mem
#BSUB -R rusage[mem=60000]
#BSUB -n 1
#BSUB -o figS1b.log
#BSUB -e figS1b.err

# FigS1B = pbmc10k KNN RNA-label-prediction accuracy vs k (reverse of FigS1A). Step 1 computes the reverse
# per-cell-type accuracy for all 23 methods (18 original latents on ${TOOLS_ROOT}/..., 5 new in the repo)
# + a forward spot-check vs published FigS1A; step 2 renders figS1b_knn_accuracy.pdf/png (same palette/order as
# FigS1A). benchmark_env has pandas/sklearn/seaborn/matplotlib. Repo-relative paths -> no edits needed.
module load conda3/202303
conda activate benchmark_env

cd ${PROJECT_ROOT}/pbmc/pbmc10k/benchmark

python 04_compute_ktest_reverse.py
python plot_figS1b.py

conda deactivate
echo "DONE: figS1b_ktest_long.csv + figS1b_knn_accuracy.pdf/png"
