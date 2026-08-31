#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J rms_figS4a_plot
#BSUB -q standard
#BSUB -R rusage[mem=16000]
#BSUB -n 1
#BSUB -o figS4a_plot.log
#BSUB -e figS4a_plot.err

# RMS Mast607A composite-metrics table (Fig5a-equivalent) from this run's staged rep1 latents.
# SPLICE_PUBLISHED=0 is REQUIRED: the splice branch in plot_metrics_matrix.R points at pbmc3k baselines
# (old/pbmc3k/...; the R script was cloned from pbmc_parse and never repointed), so splicing would inject
# pbmc3k scores into the RMS table. =0 scores all 14 RMS methods freshly. Run AFTER 00_compute_metrics_sub.sh.
module load conda3/202311
conda activate seurat4
cd ${PROJECT_ROOT}/RMS/benchmark/figS4a
export SPLICE_PUBLISHED=0
Rscript plot_metrics_matrix.R . .
conda deactivate
echo "DONE: metrics_matrix.csv + sum_metrics_clean.csv + metrics_table.html (+ .pdf/.png if webshot2 present)"
