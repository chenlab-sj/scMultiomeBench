#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J parse_umap
#BSUB -q large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o parse_umap.log
#BSUB -e parse_umap.err

# Per-method UMAP panels (one PNG per method). Run AFTER compute_metrics (needs knn_pred_label__<m>.csv).
# benchmark_env. Lists ALL 14 methods; scVI/cobolt use the normalized *_latent.csv (same as metrics).
module load conda3/202303
conda activate benchmark_env

ROOT=${PROJECT_ROOT}
cd "${ROOT}/pbmc_parse/benchmark/metrics"
LABEL=label.csv
P=${ROOT}/pbmc_parse

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
  pred="knn_pred_label__${name}.csv"
  if [ ! -f "$pred" ]; then echo "SKIP $name (run compute_metrics first: missing $pred)"; continue; fi
  echo "================ umap $name ================"
  python plot_umap.py --latent "$path" --method "$name" --label "$LABEL" --knn-pred "$pred" --out .
done

conda deactivate
echo "DONE: <method>_umap.png"
