#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J parse_metrics
#BSUB -q large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o parse_metrics.log
#BSUB -e parse_metrics.err

# Cross-platform (Parse RNA vs pbmc3k ATAC) integration metrics + KNN k=10 per method, then the
# random-adjusted ATAC accuracy. benchmark_env (scib + benchmark_fun).
# NB: this data is UNPAIRED -> 00_compute_metrics.py returns NaN for ari / ami / ks.statistic.
module load conda3/202303
conda activate benchmark_env

ROOT=${PROJECT_ROOT}
cd "${ROOT}/pbmc_parse/benchmark/metrics"
export BENCHMARK_FUN_DIR=${BENCHMARK_FUN_DIR}
LABEL=label.csv
P=${ROOT}/pbmc_parse

# name|latent_path  (standardized latent.csv; scVI/cobolt = pre-normalized by 01_prep_latents.py)
methods=(
  "BindSC|${P}/bindSC/latent.csv"
  "Seurat(CCA)|${P}/Seurat3/latent.csv"
  "scDART|${P}/scDART/latent.csv"
  "simba|${P}/simba/latent.csv"
  "scglue|${P}/scglue/latent.csv"
  "Portal|${P}/portal/latent.csv"
  "scJoint|${P}/scJoint/latent.csv"
  "scBridge|${P}/scBridge/latent.csv"
  "Cobolt|${P}/cobolt/cobolt_latent.csv"
  "scglue(multiome)|${P}/scglue_paired/latent.csv"
  "scVI|${P}/scVI/scvi_latent.csv"
  "MaxFuse|${P}/MaxFuse/latent.csv"
  "MIDAS|${P}/MIDAS/latent.csv"
  "scButterfly|${P}/scButterfly/latent.csv"
)

for entry in "${methods[@]}"; do
  name="${entry%%|*}"; path="${entry#*|}"
  if [ ! -f "$path" ]; then echo "SKIP (missing): $name  $path"; continue; fi
  echo "================ $name ================"
  python 00_compute_metrics.py --latent "$path" --method "$name" --label "$LABEL" --out .
done

conda deactivate
echo "DONE: sum_metrics.csv, celltype_metrics.csv, knn_pred_accu.csv, knn_pred_label__*.csv"
echo "NEXT: bsub < 02_adjust_sub.sh   ->  adj_atac_predaccu.csv"
