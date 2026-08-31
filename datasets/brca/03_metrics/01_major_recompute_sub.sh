#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J figS4a_major
#BSUB -q large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o figS4a_major.log
#BSUB -e figS4a_major.err

# FigS4A label-fix: recompute ONLY the 3 new methods (MaxFuse/MIDAS/scButterfly) on the MAJOR-type
# label (label_major.csv = the 7 >100-cell types, matching the published HT243 benchmark), writing into
# a separate major/ subdir so the existing 12-type CSVs are untouched. The 11 published methods are
# spliced from old/HT243B1-S1H4/ by the plot step, so they are NOT recomputed here.
# scib metrics only (sum_metrics + celltype_metrics + adj_atac_predaccu); peakdist for the new methods
# is restricted from peak/peak_similarity.csv on the laptop (no heavy peak rerun needed).
module load conda3/202303
conda activate benchmark_env

ROOT=${PROJECT_ROOT}
cd "${ROOT}/BRCA/benchmark/figS4a"
export BENCHMARK_FUN_DIR=${BENCHMARK_FUN_DIR}
LABEL=label_major.csv
mkdir -p major

methods=(
  "MaxFuse|${ROOT}/BRCA/HT243B1-S1H4/maxfuse/latent.csv"
  "MIDAS|${ROOT}/BRCA/HT243B1-S1H4/midas/latent.csv"
  "scButterfly|${ROOT}/BRCA/HT243B1-S1H4/scbutterfly/latent.csv"
)
for entry in "${methods[@]}"; do
  name="${entry%%|*}"; path="${entry#*|}"
  [ -f "$path" ] || { echo "SKIP (missing): $name $path"; continue; }
  echo "================ $name (major label) ================"
  python 00_compute_metrics.py --latent "$path" --method "$name" --label "$LABEL" --out major || echo "FAILED: $name"
done

python 01_adjust_accuracy.py --label "$LABEL" --out major   # adj_atac_predaccu for the 3 new (major-only celltypes)

conda deactivate
echo "DONE: major/{sum_metrics,celltype_metrics,adj_atac_predaccu}.csv for the 3 new methods (major-type basis)"
