#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P benchmark
#BSUB -J pbmc3k_ksens
#BSUB -q large_mem
#BSUB -R rusage[mem=32000]
#BSUB -n 1
#BSUB -o ksens.log
#BSUB -e ksens.err

# pbmc3k k-sensitivity (KNN ATAC-label-prediction accuracy vs k) with the 5 new methods, like pbmc10k FigS1A.
# Step 1 computes the 5 new methods' per-cell-type accuracy at k=5/10/20/40/80 from their repo latents; step 2
# splices them onto the published 18-method kNN_celltype_accu_sum.csv and renders pbmc3k_ksensitivity.pdf/png.
# benchmark_env has pandas/sklearn/seaborn/matplotlib. Repo-relative paths -> no edits.
module load conda3/202303
conda activate benchmark_env

cd ${PROJECT_ROOT}/pbmc/pbmc3k/benchmark/ksensitivity

python 01_compute_ktest_new_methods.py
python plot_ksens.py

conda deactivate
echo "DONE: pbmc3k_ksensitivity.pdf/png"
