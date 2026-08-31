#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J figS1c_peak
#BSUB -q large_mem
#BSUB -n 4
#BSUB -R rusage[mem=120000]
#BSUB -o figS1c_peak.log
#BSUB -e figS1c_peak.err

# FigS1C adjusted peak score (peakdist_adj). SPLICE: the 18 original methods keep the copied
# peak values (old/pbmc10k/peakdist_adj_random.csv); here we export group bigWigs + score ONLY
# the 5 NEW methods, reusing the cached annotation/random tracks in PEAK_OUT (so same baseline).
# Then combine old(18) + new(5) -> peakdist.csv (all 23) for plot_metrics_matrix.R.
# Run AFTER 00_compute_metrics_sub.sh (needs knn_pred_label__<new>.csv). Needs an R env with Signac/
# Seurat/rtracklayer/GenomicRanges/AnnotationHub/proxy/lsa/dplyr + the cached hg38 gap objects.
module load conda3/202210
conda activate seurat4

ROOT=${PROJECT_ROOT}
FIG=${ROOT}/pbmc/pbmc10k/benchmark
cd "${FIG}"

export EXPORT_BWG=${BENCHMARK_FUN_DIR}/peak_similarity/export_groupbwg.R
export PEAK_DATA=${DATA_ROOT}/pbmc10k
export FIG2B_DIR="${FIG}"
export PEAK_OUT=${BENCHMARK_FUN_DIR}/peak_similarity/pbmc10k    # has cached annotation/random + 18 old bigWigs
export NEW_METHODS=MaxFuse,MIDAS,scButterfly,MIRA,Multigrate

Rscript 02_peak_similarity.R          # -> peakdist.csv (5 new methods)

# combine copied old peak (18) + new (5) -> peakdist.csv (all 23)
Rscript -e '
suppressMessages(library(dplyr))
new <- read.csv("peakdist.csv")                                            # method, peakdist_adj (5 new)
old <- read.csv("old/pbmc10k/peakdist_adj_random.csv")                     # Method.x, peakdist_adj (18 + random)
old$method <- recode(old$Method.x, "Seurat.CCA."="Seurat(CCA)",
                     "Seurat.WNN."="Seurat(WNN)", "scglue.multiome."="scglue(multiome)")
old <- old %>% filter(Method.x != "random") %>% select(method, peakdist_adj)
allm <- bind_rows(old, new %>% select(method, peakdist_adj))
write.csv(allm, "peakdist.csv", row.names = FALSE)
cat("combined peakdist.csv:", nrow(allm), "methods (", nrow(old), "old +", nrow(new), "new )\n")
'
conda deactivate
echo "DONE: peakdist.csv (all 23 methods)"
