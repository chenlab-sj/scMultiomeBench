#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P benchmark
#BSUB -J brca_figS4b_peak
#BSUB -q large_mem
#BSUB -n 4
#BSUB -R rusage[mem=150000]
#BSUB -o figS4b_peak.log
#BSUB -e figS4b_peak.err

# Adjusted peak score (peakdist_adj) -> peakdist.csv for BRCA FigS4B (HT137). SLOW: builds
# group-coverage bigWig tracks per method x cell type from the ATAC fragments (cached -> only new
# groups re-export). seurat4 env (Signac+Seurat+rtracklayer+GenomicRanges+proxy+lsa+dplyr). Run AFTER compute.
module load conda3/202210
conda activate seurat4

ROOT=${PROJECT_ROOT}
FIG="${ROOT}/BRCA/benchmark/figS4b"
HT="${ROOT}/BRCA/HT137B1-S1H7"
cd "${FIG}"

export EXPORT_BWG=${BENCHMARK_FUN_DIR}/peak_similarity/export_groupbwg.R
export PEAK_DATA="${HT}"
export FRAG_FILE="${HT}/HT137B1-S1H7-fragments.tsv.gz"   # 1.6G fragments (+ .tbi)
export H5_FILE="${HT}/HT137B1-S1H7_commonpeak"           # test ATAC, common peaks MTX dir (the .h5 was written corrupt; the dir is complete)
export FIG2B_DIR="${FIG}"                                # holds knn_pred_label__*.csv + label.csv
export PEAK_OUT="${FIG}/peak"

Rscript peak_similarity.R

conda deactivate
echo "DONE: peakdist.csv (method, peakdist_adj)"
