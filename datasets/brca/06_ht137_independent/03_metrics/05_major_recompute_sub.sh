#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

#BSUB -P benchmark
#BSUB -J figS4b_major
#BSUB -q large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o figS4b_major.log
#BSUB -e figS4b_major.err

# FigS4B (HT137) label-fix. This is a NEW (revision) analysis with NO published values to splice, so
# ALL 14 methods are recomputed on label_major.csv (the >100-cell types) into major/, giving a
# major-type-consistent composite (matching figS4a's corrected basis for the R1 comparison).
module load conda3/202303
conda activate benchmark_env

ROOT=${PROJECT_ROOT}
cd "${ROOT}/BRCA/benchmark/figS4b"
export BENCHMARK_FUN_DIR=${BENCHMARK_FUN_DIR}
LABEL=label_major.csv
HT="${ROOT}/BRCA/HT137B1-S1H7"
mkdir -p major

methods=(
  "Seurat(CCA)|${HT}/seurat3/coembed_coor.csv"
  "BindSC|${HT}/bindsc/coembed_coor.csv"
  "scBridge|${HT}/scBridge/latent.csv"
  "scglue|${HT}/scglue/latent.csv"
  "scglue(multiome)|${HT}/scglue_paired/latent.csv"
  "scJoint|${HT}/scjoint/latent.csv"
  "scVI|${HT}/scVI/rep2/scvi_latent.csv"
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
  [ -f "$path" ] || { echo "SKIP (missing): $name $path"; continue; }
  echo "================ $name (major label) ================"
  python compute_metrics.py --latent "$path" --method "$name" --label "$LABEL" --out major || echo "FAILED: $name"
done

python adjust_accuracy.py --label "$LABEL" --out major

conda deactivate
echo "DONE: major/{sum_metrics,celltype_metrics,adj_atac_predaccu}.csv for all 14 methods (major-type basis)"
