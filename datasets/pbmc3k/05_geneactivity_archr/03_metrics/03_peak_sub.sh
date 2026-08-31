#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P benchmark
#BSUB -J ga_archr_peak
#BSUB -q large_mem
#BSUB -n 4
#BSUB -R rusage[mem=150000]
#BSUB -o ga_archr_peak.log
#BSUB -e ga_archr_peak.err
# pbmc3k peak score. 02_peak_similarity.R defaults FRAG/H5 to the pbmc3k filenames under PEAK_DATA.
module load conda3/202210
conda activate seurat4
cd ${PROJECT_ROOT}/pbmc/pbmc3k/benchmark/geneactivity_archr/metrics
export EXPORT_BWG=${BENCHMARK_FUN_DIR}/peak_similarity/export_groupbwg.R
export PEAK_DATA=${DATA_ROOT}/pbmc3k
export FIG2B_DIR=${PROJECT_ROOT}/pbmc/pbmc3k/benchmark/geneactivity_archr/metrics
export PEAK_OUT=${PROJECT_ROOT}/pbmc/pbmc3k/benchmark/geneactivity_archr/metrics/peak
Rscript 02_peak_similarity.R
conda deactivate
echo "DONE: peakdist.csv"
