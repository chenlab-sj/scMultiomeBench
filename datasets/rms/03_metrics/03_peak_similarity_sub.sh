#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J rms_figS4a_peak
#BSUB -q large_mem
#BSUB -n 4
#BSUB -R rusage[mem=150000]
#BSUB -o figS4a_peak.log
#BSUB -e figS4a_peak.err

# RMS Mast607A peak similarity (peakdist_adj) for Fig5a. Recomputes ALL 14 methods consistently from this
# run's knn_pred_label__*.csv (written by 00_compute_metrics_sub.sh), using the hg38 Mast607A ATAC -- matching
# the original peak setup (genome hg38, chr-prefixed peaks, hg38 gap tracks). The bigWig export is the slow
# part and is CACHED (a group whose .bw exist is skipped). Run AFTER 00_compute_metrics_sub.sh, BEFORE plot_metrics_matrix_sub.sh.
module load conda3/202210
conda activate seurat4

ROOT=${PROJECT_ROOT}
cd "${ROOT}/RMS/benchmark/figS4a"
DATA38="${ROOT}/RMS/Mast607/data/Mast607A_TB19_22652/hg38/outs"
export EXPORT_BWG=${BENCHMARK_FUN_DIR}/export_groupbwg.R
export FRAG="${DATA38}/atac_fragments.tsv.gz"
export H5="${DATA38}/filtered_feature_bc_matrix.h5"
export FIG2B_DIR="${ROOT}/RMS/benchmark/figS4a"     # holds knn_pred_label__*.csv + label.csv
export PEAK_OUT="${ROOT}/RMS/benchmark/figS4a/peak"

Rscript 02_peak_similarity.R

conda deactivate
echo "DONE: peakdist.csv (14 methods incl. MaxFuse/MIDAS/scButterfly)"
