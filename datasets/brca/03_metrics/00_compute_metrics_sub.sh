#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J brca_figS4a_metrics
#BSUB -q large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o figS4a_metrics.log
#BSUB -e figS4a_metrics.err

# BRCA HT243B1-S1H4 benchmark matrix (FigS4A). Same parameterized 00_compute_metrics.py +
# 01_adjust_accuracy.py as pbmc3k fig2b/, just BRCA's label.csv + method->latent list.
# 18 existing methods read their cluster latents (${DATA_ROOT}/HTAN/HT243B1-S1H4/<m>/);
# 3 new methods read the repo latents. Missing/format-bad latents are SKIPPED (loop continues) ->
# iterate per method like pbmc3k did for Conos.
module load conda3/202303
conda activate benchmark_env

ROOT=${PROJECT_ROOT}
cd "${ROOT}/BRCA/benchmark/figS4a"
export BENCHMARK_FUN_DIR=${BENCHMARK_FUN_DIR}
LABEL=label.csv
HTAN=${DATA_ROOT}/HTAN/HT243B1-S1H4

# FigS4A = TOP-PERFORMER subset (14 methods): the 11 carried into the deep-dive figures
# (the set used in HT243_S1H4_reproduce.ipynb) + the 3 new methods. EXCLUDED below-cutoff:
# Conos, LIGER, MinNet, scMoMaT, Unioncom, MultiMAP, Seurat(WNN). Paths from benchmark.ipynb.
methods=(
  "Seurat(CCA)|${HTAN}/Seurat3/coembed_coor.csv"
  "BindSC|${HTAN}/bindSC/coembed_coor.csv"
  "scBridge|${HTAN}/scBridge/latent.csv"
  "scglue|${HTAN}/scglue/latent.csv"
  "scglue(multiome)|${HTAN}/scglue_paired/scglue_latent.csv"
  "scJoint|${HTAN}/scJoint/latent.csv"
  "scVI|${HTAN}/scVI/scvi_latent.csv"                                  # CORRECT path (the figure scripts had HT514B1-S1H3 = wrong patient)
  "Cobolt|${HTAN}/cobolt/cobolt_latent.csv"
  "simba|${HTAN}/simba/latent.csv"
  "Portal|${HTAN}/portal/lat_df.csv"
  "scDART|${HTAN}/scDART/latent.csv"
  # ---- new methods (this revision) ----
  "MaxFuse|${ROOT}/BRCA/HT243B1-S1H4/maxfuse/latent.csv"
  "MIDAS|${ROOT}/BRCA/HT243B1-S1H4/midas/latent.csv"
  "scButterfly|${ROOT}/BRCA/HT243B1-S1H4/scbutterfly/latent.csv"
)

for entry in "${methods[@]}"; do
  name="${entry%%|*}"; path="${entry#*|}"
  if [ ! -f "$path" ]; then echo "SKIP (missing): $name  $path"; continue; fi
  echo "================ $name ================"
  python 00_compute_metrics.py --latent "$path" --method "$name" --label "$LABEL" --out . || echo "FAILED: $name (iterate on its format)"
done

python 01_adjust_accuracy.py --label "$LABEL" --out .

conda deactivate
echo "DONE: sum_metrics.csv, celltype_metrics.csv, knn_pred_accu.csv, adj_atac_predaccu.csv"
