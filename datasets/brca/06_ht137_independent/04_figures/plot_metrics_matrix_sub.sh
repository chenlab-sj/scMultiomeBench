#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P benchmark
#BSUB -J brca_figS4b_plot
#BSUB -q large_mem
#BSUB -R rusage[mem=40000]
#BSUB -n 1
#BSUB -o figS4b_plot.log
#BSUB -e figS4b_plot.err

# BRCA FigS4B matrix (= original figS4c, HT137, GENUINELY UNPAIRED -> plot_metrics.R drops the same-cell
# columns ks.statistic/ARI/AMI). SPLICE_PUBLISHED=1 -> the 11 OLD methods keep their ORIGINAL published
# figS4c values (old/HT137B1-S1H7/ + old peak peakdist_adj_random.csv); only the 3 NEW methods
# (MIDAS/MaxFuse/scButterfly, count>100 major basis from major/) are added, then re-ranked.
# REQUIRES: old/HT137B1-S1H7/peakdist_adj_random.csv (old HT137 peak; Method.x + peakdist_adj).
# Needs R with dplyr + formattable (+ webshot2 for the PNG).
module load conda3/202311
conda activate seurat4

ROOT=${PROJECT_ROOT}
cd "${ROOT}/BRCA/benchmark/figS4b"

export SPLICE_PUBLISHED=1
Rscript plot_metrics.R . .

conda deactivate
echo "DONE: metrics_matrix.csv (= BRCA FigS4B) + sum_metrics_clean.csv + metrics_table.html/png"
