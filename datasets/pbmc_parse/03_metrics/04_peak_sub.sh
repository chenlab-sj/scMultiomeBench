#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J parse_peak
#BSUB -q large_mem
#BSUB -n 4
#BSUB -R rusage[mem=120000]
#BSUB -o parse_peak.log
#BSUB -e parse_peak.err

# Adjusted peak score (peakdist_adj) -> peakdist.csv. SLOW first run: builds group-coverage bigWig
# tracks per method x cell type (cached -> only new groups re-export). The ATAC is the SAME pbmc3k
# ATAC used in metrics, so PEAK_DATA is unchanged. Run AFTER compute_metrics (needs knn_pred_label__*.csv).
# Needs an R env with Signac + Seurat + rtracklayer + GenomicRanges + AnnotationHub + proxy + lsa + dplyr.
module load conda3/202210
conda activate seurat4          # <-- adjust to an R env with the Bioconductor deps above

ROOT=${PROJECT_ROOT}
METRICS="${ROOT}/pbmc_parse/benchmark/metrics"
cd "${METRICS}"

export EXPORT_BWG=${BENCHMARK_FUN_DIR}/peak_similarity/export_groupbwg.R
export PEAK_DATA=${DATA_ROOT}/pbmc3k        # same pbmc3k ATAC fragments .tsv.gz + test .h5
export METRICS_DIR="${METRICS}"
export PEAK_OUT="${METRICS}/peak"

Rscript 02_peak_similarity.R

conda deactivate
echo "DONE: peakdist.csv (method, peakdist_adj)"
