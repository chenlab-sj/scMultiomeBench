#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J brca_fig4_nb
#BSUB -q large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o fig4_reproduce_nb.log
#BSUB -e fig4_reproduce_nb.err

# Regenerate BRCA Fig4 reproducibility on the 7-type label (11 old + 3 new methods), manuscript style.
# Two steps: (1) run 03_fig4_reproduce_7type.ipynb to compute the metrics, (2) render the manuscript-style panels.
# Manuscript panels:
#   Fig4A = pairwise louvain NMI boxplot   -> fig4A_nmi_boxplot_7type.pdf/png
#   Fig4B = per-cell-type SD dot+bar       -> fig4B_sd_dotbar_7type.pdf/png
# Re-selected reps: MIDAS = base/rep3/rep4, Cobolt = base/rep2/rep5 (from midas_pick.csv / cobolt_pick.csv);
# the other 12 methods are unchanged. The 10 UNCHANGED old methods reproduce the PUBLISHED Fig4A values
# exactly; Cobolt is deliberately re-selected so its value changes. The SD panel is CPU-only; the NMI panel
# needs louvain clustering of the 4 re-selected/new methods (MaxFuse, MIDAS, scButterfly, Cobolt: per-rep
# find_louvain_res -> nclust=7, seed 420; scanpy + benchmark_fun, hence this job). The other 10 old NMI
# values reuse the published 7-type clusters.
module load conda3/202303
conda activate benchmark_env

ROOT=${PROJECT_ROOT}
cd "${ROOT}/BRCA/benchmark/fig4"
export BENCHMARK_FUN_DIR=${BENCHMARK_FUN_DIR}

# 1) run the notebook -> reproducibility_7type.csv, nmi_df_7type.csv, rep_knn_k10_pred_label_7type.csv
#    (Fig4B new-method louvain clustering happens here, so nmi_df now has all 14 methods)
jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=3600 03_fig4_reproduce_7type.ipynb

# 2) render the manuscript-style panels from those CSVs (all 14 methods in BOTH panels now):
#    fig4A_nmi_boxplot_7type.pdf/png (NMI boxplot) + fig4B_sd_dotbar_7type.pdf/png (per-cell-type SD dot+bar)
python plot_fig4_manuscript_7type.py

conda deactivate
echo "DONE: Fig4A (NMI) + Fig4B (SD dot+bar), manuscript style, all 14 methods incl. MaxFuse/MIDAS/scButterfly"
