#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J bmmc_xs_major
#BSUB -q large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o cs_major.log
#BSUB -e cs_major.err

# BMMC cross-site (step 2, after compute). NEW analysis, no published values to splice -> ALL 12 methods
# recomputed on label_major.csv (test-ATAC > 100-cell types) into major/, then adjust on the major basis.
# plot_metrics_matrix.R (SPLICE_PUBLISHED=0) reads major/{sum,celltype,adj_atac_predaccu}.csv for the composite
# score, plus the top-level peakdist.csv. Uses the SAME normalized latents/<key>.csv prep_latents wrote.
# LS_SUBCWD = the experiment dir -> one script serves both experiments.
module load conda3/202303
conda activate benchmark_env

FIG="${LS_SUBCWD:-$PWD}"
cd "${FIG}" || exit 1
export BENCHMARK_FUN_DIR=${BENCHMARK_FUN_DIR}
LABEL=label_major.csv
mkdir -p major

methods=(
  "Seurat(CCA)|latents/seurat3.csv"
  "BindSC|latents/bindsc.csv"
  "scBridge|latents/scBridge.csv"
  "scglue|latents/scglue.csv"
  "scglue(multiome)|latents/scglue_paired.csv"
  "scJoint|latents/scjoint.csv"
  "scVI|latents/scVI.csv"
  "simba|latents/simba.csv"
  "Portal|latents/portal.csv"
  "MaxFuse|latents/maxfuse.csv"
  "MIDAS|latents/midas.csv"
  "scButterfly|latents/scbutterfly.csv"
)
for entry in "${methods[@]}"; do
  name="${entry%%|*}"; path="${entry#*|}"
  [ -f "$path" ] || { echo "SKIP (missing): $name $path"; continue; }
  echo "================ $name (major label) ================"
  python 00_compute_metrics.py --latent "$path" --method "$name" --label "$LABEL" --out major || echo "FAILED: $name"
done

python 01_adjust_accuracy.py --label "$LABEL" --out major

conda deactivate
echo "DONE: major/{sum_metrics,celltype_metrics,adj_atac_predaccu}.csv (major-type basis)"
