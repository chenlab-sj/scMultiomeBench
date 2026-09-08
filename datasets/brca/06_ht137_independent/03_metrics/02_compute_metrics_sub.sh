#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P benchmark
#BSUB -J brca_figS4b_metrics
#BSUB -q large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o figS4b_metrics.log
#BSUB -e figS4b_metrics.err

# BRCA HT137B1-S1H7 benchmark matrix (Fig S4B) -- GENUINELY UNPAIRED, so this mirrors the Parse
# cross-platform pipeline: same parameterized compute_metrics.py, just HT137's label.csv + latents.
# Step 0 (re)builds label.csv and normalizes the scVI latents, then scores all 14 methods.
# Missing/format-bad latents are SKIPPED (loop continues) -- e.g. a new method whose run hasn't
# finished yet. The 2 scVI seed rows are DIAGNOSTIC (your "try another seed"): they let you compare
# scVI at seeds 420/0/40 -> pick the representative one for the headline 14-method figure.
module load conda3/202303
conda activate benchmark_env

ROOT=${PROJECT_ROOT}
FIG="${ROOT}/BRCA/benchmark/figS4b"
HT="${ROOT}/BRCA/HT137B1-S1H7"
cd "${FIG}"
export BENCHMARK_FUN_DIR=${BENCHMARK_FUN_DIR}
LABEL=label.csv

# ---- step 0: build label.csv (harmonize ATAC->RNA vocab) + normalize scVI latents ----
python build_label.py
python prep_latents.py

# ---- the 14 benchmark methods (same set as HT243 figS4a) ----
methods=(
  "Seurat(CCA)|${HT}/seurat3/coembed_coor.csv"
  "BindSC|${HT}/bindsc/coembed_coor.csv"
  "scBridge|${HT}/scBridge/latent.csv"
  "scglue|${HT}/scglue/latent.csv"
  "scglue(multiome)|${HT}/scglue_paired/latent.csv"
  "scJoint|${HT}/scjoint/latent.csv"
  "scVI|${HT}/scVI/rep2/scvi_latent.csv"        # seed 0 = best of the 3 scVI seeds (0.67 vs 0.58/0.59)
  "Cobolt|${HT}/cobolt/cobolt_latent.csv"
  "simba|${HT}/simba/latent.csv"
  "Portal|${HT}/portal/lat_df.csv"
  "scDART|${HT}/scDART/latent.csv"
  # ---- new methods (this revision) ----
  "MaxFuse|${HT}/maxfuse/latent.csv"
  "MIDAS|${HT}/midas/latent.csv"
  "scButterfly|${HT}/scbutterfly/latent.csv"
)

for entry in "${methods[@]}"; do
  name="${entry%%|*}"; path="${entry#*|}"
  if [ ! -f "$path" ]; then echo "SKIP (missing): $name  $path"; continue; fi
  echo "================ $name ================"
  python compute_metrics.py --latent "$path" --method "$name" --label "$LABEL" --out . || echo "FAILED: $name (iterate on its format)"
done

conda deactivate
echo "DONE: sum_metrics.csv, celltype_metrics.csv, knn_pred_label__*.csv (then: adjust + peak -> plot)"
