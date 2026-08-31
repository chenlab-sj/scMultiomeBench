#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J fig2b_umap
#BSUB -q large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o fig2b_umap.log
#BSUB -e fig2b_umap.err

# Per-method UMAP panels. Run AFTER 00_compute_metrics_sub.sh (needs knn_pred_label__<m>.csv).
# Edit the list to just the methods you want figures for.
module load conda3/202303
conda activate benchmark_env

ROOT=${PROJECT_ROOT}
cd "${ROOT}/pbmc/pbmc3k/benchmark/fig2b"
LABEL=label.csv

methods=(
  "MaxFuse|${ROOT}/pbmc/pbmc3k/scripts/maxfuse/latent.csv"
  "Multigrate|${ROOT}/pbmc/pbmc3k/scripts/multigrate/latent.csv"
  "MIDAS|${ROOT}/pbmc/pbmc3k/scripts/midas/latent.csv"
  "scButterfly|${ROOT}/pbmc/pbmc3k/scripts/scbutterfly/latent.csv"
  "MIRA|${ROOT}/pbmc/pbmc3k/scripts/MIRA/latent.csv"
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
