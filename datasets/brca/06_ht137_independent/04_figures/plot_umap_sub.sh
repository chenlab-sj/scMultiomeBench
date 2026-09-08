#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P benchmark
#BSUB -J brca_figS4b_umap
#BSUB -q large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o figS4b_umap.log
#BSUB -e figS4b_umap.err

# Per-method UMAP panels (one PNG per method). Run AFTER compute (needs knn_pred_label__<m>.csv).
# benchmark_env. Same 14 methods as the matrix; scVI/cobolt use the normalized *_latent.csv.
module load conda3/202303
conda activate benchmark_env

ROOT=${PROJECT_ROOT}
FIG="${ROOT}/BRCA/benchmark/figS4b"
HT="${ROOT}/BRCA/HT137B1-S1H7"
cd "${FIG}"
LABEL=label.csv

methods=(
  "Seurat(CCA)|${HT}/seurat3/coembed_coor.csv"
  "BindSC|${HT}/bindsc/coembed_coor.csv"
  "scBridge|${HT}/scBridge/latent.csv"
  "scglue|${HT}/scglue/latent.csv"
  "scglue(multiome)|${HT}/scglue_paired/latent.csv"
  "scJoint|${HT}/scjoint/latent.csv"
  "scVI|${HT}/scVI/scvi_latent.csv"
  "Cobolt|${HT}/cobolt/cobolt_latent.csv"
  "simba|${HT}/simba/latent.csv"
  "Portal|${HT}/portal/lat_df.csv"
  "scDART|${HT}/scDART/latent.csv"
  "MaxFuse|${HT}/maxfuse/latent.csv"
  "MIDAS|${HT}/midas/latent.csv"
  "scButterfly|${HT}/scbutterfly/latent.csv"
)

for entry in "${methods[@]}"; do
  name="${entry%%|*}"; path="${entry#*|}"
  pred="knn_pred_label__${name}.csv"
  if [ ! -f "$path" ]; then echo "SKIP $name (missing latent: $path)"; continue; fi
  if [ ! -f "$pred" ]; then echo "SKIP $name (run compute first: missing $pred)"; continue; fi
  echo "================ umap $name ================"
  python plot_umap.py --latent "$path" --method "$name" --label "$LABEL" --knn-pred "$pred" --out .
done

conda deactivate
echo "DONE: <method>_umap.png"
