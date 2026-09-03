#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J fig2b_peak
#BSUB -q large_mem
#BSUB -n 4
#BSUB -R rusage[mem=120000]
#BSUB -o fig2b_peak.log
#BSUB -e fig2b_peak.err

# Adjusted peak score (peakdist_adj) -> peakdist.csv. SLOW: builds group-coverage bigWig tracks
# per method x cell type (cached across runs -> only new groups re-export). Needs R with
# Signac + Seurat + rtracklayer + GenomicRanges + AnnotationHub + proxy + lsa + dplyr, and the
# AnnotationHub hg38 gap objects (AH107355-9) -- already cached from your original peak run.
module load conda3/202210
conda activate seurat4          # <-- adjust to an R env that has the Bioconductor deps above

ROOT=${PROJECT_ROOT}
FIG2B="${ROOT}/pbmc/pbmc3k/benchmark/fig2b"
cd "${FIG2B}"

export EXPORT_BWG=${BENCHMARK_FUN_DIR}/export_groupbwg.R
export PEAK_DATA=${DATA_ROOT}/pbmc3k        # fragments .tsv.gz + test .h5 live here
export FIG2B_DIR="${FIG2B}"
export PEAK_OUT="${FIG2B}/peak"

Rscript 02_peak_similarity.R

conda deactivate
echo "DONE: peakdist.csv (method, peakdist_adj)"
