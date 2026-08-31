#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J ga_archr_metrics
#BSUB -q large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o ga_archr_metrics.log
#BSUB -e ga_archr_metrics.err

# R2.1: score the 6 ArchR-gene-activity variant latents with the SAME paired composite as fig2b
# (pbmc3k is multiome-split = same cells, so the full composite incl. same-cell metrics is valid).
# Uses fig2b's label.csv (copied here). compute SKIPs any latent whose run hasn't finished.
module load conda3/202303
conda activate benchmark_env

GA=${PROJECT_ROOT}/pbmc/pbmc3k/benchmark/geneactivity_archr
cd "${GA}/metrics"
export BENCHMARK_FUN_DIR=${BENCHMARK_FUN_DIR}
LABEL=label.csv

methods=(
  "Seurat(CCA)|${GA}/Seurat_CCA/res_pbmc3k_testall/lat_df.csv"
  "BindSC|${GA}/bindsc/res_pbmc3k/lat_df.csv"
  "MaxFuse|${GA}/maxfuse/latent.csv"
  "scBridge|${GA}/scBridge/latent.csv"
  "Portal|${GA}/portal/lat_df.csv"
  "scJoint|${GA}/scJoint/output/scJoint_latent.csv"
)

for entry in "${methods[@]}"; do
  name="${entry%%|*}"; path="${entry#*|}"
  if [ ! -f "$path" ]; then echo "SKIP (missing): $name  $path"; continue; fi
  echo "================ $name ================"
  python 00_compute_metrics.py --latent "$path" --method "$name" --label "$LABEL" --out . || echo "FAILED: $name"
done

conda deactivate
echo "DONE: sum_metrics.csv, celltype_metrics.csv, knn_pred_label__*.csv (then adjust + peak -> plot)"
