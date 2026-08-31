#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J brca_figS4a_peak
#BSUB -q large_mem
#BSUB -n 4
#BSUB -R rusage[mem=150000]
#BSUB -o figS4a_peak.log
#BSUB -e figS4a_peak.err

# Adjusted peak score (peakdist_adj) -> peakdist.csv for BRCA FigS4A. SLOW: builds group-coverage
# bigWig tracks per method x cell type from the ATAC fragments (cached -> only new groups re-export).
# Needs R with Signac+Seurat+rtracklayer+GenomicRanges+AnnotationHub+proxy+lsa+dplyr and the cached
# hg38 gap objects (already cached from your original peak run).
module load conda3/202210
conda activate seurat4

ROOT=${PROJECT_ROOT}
FIG="${ROOT}/BRCA/benchmark/figS4a"
BRCA="${ROOT}/BRCA/HT243B1-S1H4"
cd "${FIG}"

export EXPORT_BWG=${BENCHMARK_FUN_DIR}/peak_similarity/export_groupbwg.R
export PEAK_DATA="${BRCA}"
export FRAG_FILE="${BRCA}/HT243B1-S1H4-atac_fragments.tsv.gz"   # 2.7G fragments (+ .tbi)
export H5_FILE="${BRCA}/HT243B1-S1H4_commonpeaks.h5"            # test ATAC, common peaks (matches fragments)
export FIG2B_DIR="${FIG}"                                       # holds knn_pred_label__*.csv + label.csv
export PEAK_OUT="${FIG}/peak"

Rscript 02_peak_similarity.R

conda deactivate
echo "DONE: peakdist.csv (method, peakdist_adj)"
