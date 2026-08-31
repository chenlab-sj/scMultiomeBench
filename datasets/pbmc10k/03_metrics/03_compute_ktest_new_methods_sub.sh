#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P benchmark
#BSUB -J figS1a
#BSUB -q large_mem
#BSUB -R rusage[mem=32000]
#BSUB -n 1
#BSUB -o figS1a.log
#BSUB -e figS1a.err

# FigS1A = pbmc10k KNN ATAC-label-prediction accuracy vs k value, all methods (18 old + 5 new).
# Step 1 computes the 5 new methods' per-k per-celltype accuracy from their latents (skips any not run yet);
# step 2 renders figS1a_knn_accuracy.pdf/png (auto-includes whatever methods are available).
# benchmark_env has pandas/sklearn/seaborn/matplotlib. Scripts use repo-relative paths -> no path edits needed.
module load conda3/202303
conda activate benchmark_env

cd ${PROJECT_ROOT}/pbmc/pbmc10k/benchmark

python 01_compute_ktest_new_methods.py
python plot_figS1a.py

conda deactivate
echo "DONE: figS1a_knn_accuracy.pdf/png"
