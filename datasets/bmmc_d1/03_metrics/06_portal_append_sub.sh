#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P benchmark
#BSUB -J bmmc_portal
#BSUB -q large_mem
#BSUB -R rusage[mem=200000]
#BSUB -n 4
#BSUB -o portal.log
#BSUB -e portal.err
# Re-add Portal to Fig6 by computing ONLY Portal's metrics and APPENDING them to the existing CSVs --
# the 18 already-scored methods are NOT re-run. Prereqs already done on the Mac side:
#   staged/Portal/latent.csv, knn_pred_label__Portal.csv, and Portal's rows in knn_pred_accu.csv +
#   adj_atac_predaccu.csv. This job fills the two cluster-only pieces: the scib metrics and the peak step.
set -e
ROOT=${PROJECT_ROOT}
cd "${ROOT}/BMMC_d1/benchmark"
export BENCHMARK_FUN_DIR=${BENCHMARK_FUN_DIR}

# 1) scib batch-aware metrics -> APPEND Portal row to sum_metrics.csv + celltype_metrics.csv
#    (PORTAL_ONLY=1 scores just Portal into portal_only/ then appends; ~1/18 of the full run's time)
module load conda3/202303
conda activate benchmark_env
PORTAL_ONLY=1 python 02_run_bmmc_metrics.py
conda deactivate

# 2) peak similarity -> APPEND Portal to peakdist.csv. PEAK_LAST=Portal keeps the existing methods' M1..M18
#    tokens (and their cached peak/M*/ bigWigs) valid, so only Portal's coverage is exported fresh.
module load conda3/202210
conda activate seurat4
export EXPORT_BWG=${BENCHMARK_FUN_DIR}/export_groupbwg.R
export FIG2B_DIR=.
export PEAK_OUT=./peak
export LABEL=./old/BMMC_d1/label.csv
export BMMC_SRA=${DATA_ROOT}/BMMC/NCBI_sra
export PEAK_LAST=Portal
Rscript 02_peak_similarity.R
conda deactivate

echo "DONE: Portal appended to sum_metrics.csv, celltype_metrics.csv, peakdist.csv."
echo "Next: re-run  python fig6/plot_fig6.py  -> Portal switches from SPLICED(published) to computed."
