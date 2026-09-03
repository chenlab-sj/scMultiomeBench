#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J bmmc_sp_peak
#BSUB -q large_mem
#BSUB -n 4
#BSUB -R rusage[mem=150000]
#BSUB -o sp_peak.log
#BSUB -e sp_peak.err

# Adjusted peak score (peakdist_adj) -> peakdist.csv for the BMMC s1d1 paired-multiome baseline. Test
# ATAC = s1d1, so this uses s1d1's fragments + the s1d1 common-peak test h5. seurat4 env, genome hg38.
# Runs after compute (reads its knn_pred_label__*.csv + label.csv). LS_SUBCWD = the experiment dir.
module load conda3/202210
conda activate seurat4

FIG="${LS_SUBCWD:-$PWD}"
cd "${FIG}" || exit 1

export EXPORT_BWG=${BENCHMARK_FUN_DIR}/export_groupbwg.R
export FRAG_FILE=${DATA_ROOT}/BMMC/NCBI_sra/s1d1/outs/atac_fragments.tsv.gz
export H5_FILE=${DATA_ROOT}/BMMC_d1/test3_s1d1.h5        # s1d1 test h5 (Gene Expression + Peaks)
export FIG2B_DIR="${FIG}"                                          # holds knn_pred_label__*.csv + label.csv
export PEAK_OUT="${FIG}/peak"

Rscript 02_peak_similarity.R

conda deactivate
echo "DONE: peakdist.csv (method, peakdist_adj)"
