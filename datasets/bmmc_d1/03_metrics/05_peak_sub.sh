#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P benchmark
#BSUB -J bmmc_peak
#BSUB -q large_mem
#BSUB -R rusage[mem=250000]
#BSUB -n 4
#BSUB -o peak.log
#BSUB -e peak.err
# BMMC Fig2b peak step -> peakdist.csv. Builds per-cell-type coverage bigWigs from the 3 test batches'
# cellranger-arc fragments (s2d1/s4d1/s1d1, ~6.4 GB total), grouped by each method's KNN-predicted
# labels + truth + random, then scaled-euclidean similarity vs truth -> peakdist_adj. HEAVY + long
# (18 methods x ~18 cell types x 3 batches of coverage export). CACHED + RESUMABLE: bigWigs already
# in peak/<group>/ are skipped, so if this dies you can just resubmit and it continues.
# Prereq: 03_knn_bmmc.py already wrote knn_pred_label__*.csv here (done).
module load conda3/202210
conda activate seurat4

ROOT=${PROJECT_ROOT}
cd "${ROOT}/BMMC_d1/benchmark"
export EXPORT_BWG=${BENCHMARK_FUN_DIR}/peak_similarity/export_groupbwg.R
export FIG2B_DIR=.
export PEAK_OUT=./peak
export LABEL=./old/BMMC_d1/label.csv
export BMMC_SRA=${DATA_ROOT}/BMMC/NCBI_sra

Rscript 02_peak_similarity.R

conda deactivate
echo "DONE: peakdist.csv (+ peak/peak_similarity.csv)"
