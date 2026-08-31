#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J brca_figS4a_plot
#BSUB -q large_mem
#BSUB -R rusage[mem=40000]
#BSUB -n 1
#BSUB -o figS4a_plot.log
#BSUB -e figS4a_plot.err

# BRCA FigS4A matrix. SPLICE_PUBLISHED=1 -> the 12 OLD methods (incl MinNet) keep their ORIGINAL
# published figS4a values (old/HT243B1-S1H4_adj/ + old peak peakdist_adj_random.csv), and only the 3
# NEW methods (MIDAS/MaxFuse/scButterfly, count>100 major basis from major/) are added, then re-ranked.
# REQUIRES: old/HT243B1-S1H4_adj/peakdist_adj_random.csv (old HT243 peak; Method.x + peakdist_adj).
# Needs R with dplyr + formattable. Writes fig2b_matrix.csv + sum_metrics_clean.csv + fig2b_table.html.
module load conda3/202311
conda activate seurat4

ROOT=${PROJECT_ROOT}
cd "${ROOT}/BRCA/benchmark/figS4a"

export SPLICE_PUBLISHED=1
Rscript plot_metrics_matrix.R . .

conda deactivate
echo "DONE: fig2b_matrix.csv (= BRCA FigS4A) + fig2b_table.html"
