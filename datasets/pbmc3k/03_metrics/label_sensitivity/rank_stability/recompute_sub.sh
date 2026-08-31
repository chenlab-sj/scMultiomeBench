#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J labelsens_recompute
#BSUB -q large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o recompute.log
#BSUB -e recompute.err

# R4 label-circularity rank-stability: re-score the SAME pbmc3k method latents (NO method re-runs) against
# the independent SingleR RNA-only labels (label_singleR.csv) instead of the 10X labels. 00_compute_metrics.py
# + 01_adjust_accuracy.py are byte-identical to fig2b's; only --label changes. Output -> this dir; then
# rank_stability.py builds the SingleR composite (reusing fig2b's peakdist) and Spearmans vs the 10X-label ranking.
module load conda3/202303
conda activate benchmark_env

ROOT=${PROJECT_ROOT}
cd "${ROOT}/pbmc/pbmc3k/benchmark/label_sensitivity/rank_stability"
export BENCHMARK_FUN_DIR=${BENCHMARK_FUN_DIR}
LABEL="${ROOT}/pbmc/pbmc3k/benchmark/label_sensitivity/label_singleR.csv"

methods=(
  "BindSC|${TOOLS_ROOT}/bindsc/pbmc3k/res_pbmc3k/coembed_coor.csv"
  "Conos|${TOOLS_ROOT}/Concos/pbmc3k/res_pbmc3k/coembed_coor.csv"
  "LIGER|${TOOLS_ROOT}/liger/pbmc3k/res_pbmc3k/coembed_coor.csv"
  "MultiMAP|${TOOLS_ROOT}/multimap/pbmc3k/res_pbmc3k/coembed_coor.csv"
  "Seurat(CCA)|${TOOLS_ROOT}/Seuratv3/pbmc3k/res_pbmc3k_testall/lat_df.csv"
  "scDART|${TOOLS_ROOT}/scDART/scDART/pbmc3k/res_pbmc3k/latent.csv"
  "scMoMaT|${TOOLS_ROOT}/scMoMaT/pbmc3k/latent.csv"
  "simba|${TOOLS_ROOT}/simba/pbmc3k/latent.csv"
  "scglue|${TOOLS_ROOT}/scglue/pbmc3k/res_pbmc3k/latent.csv"
  "Unioncom|${TOOLS_ROOT}/unioncom/pbmc3k/latent.csv"
  "Portal|${TOOLS_ROOT}/portal/Portal/pbmc3k/lat_df.csv"
  "scJoint|${TOOLS_ROOT}/scJoint/pbmc3k/output/scJoint_latent.csv"
  "scBridge|${TOOLS_ROOT}/scBridge/pbmc3k/latent.csv"
  "Cobolt|${TOOLS_ROOT}/cobolt/pbmc3k/res_pbmc3k/cobolt_latent.csv"
  "scglue(multiome)|${TOOLS_ROOT}/scglue/scglue_withpair/pbmc3k/scglue_latent.csv"
  "scVI|${TOOLS_ROOT}/scvi/pbmc3k/res_pbmc3k/scvi_latent.csv"
  "Seurat(WNN)|${TOOLS_ROOT}/Seuratv4/pbmc3k/res_pbmc3k_seurat4/seurat4_latent.csv"
  "MinNet|${TOOLS_ROOT}/MinNet/pbmc3k/coembed_coor.csv"
  "MaxFuse|${ROOT}/pbmc/pbmc3k/scripts/maxfuse/latent.csv"
  "Multigrate|${ROOT}/pbmc/pbmc3k/scripts/multigrate/latent.csv"
  "MIDAS|${ROOT}/pbmc/pbmc3k/scripts/midas/latent.csv"
  "scButterfly|${ROOT}/pbmc/pbmc3k/scripts/scbutterfly/latent.csv"
  "MIRA|${ROOT}/pbmc/pbmc3k/scripts/MIRA/latent.csv"
)

for entry in "${methods[@]}"; do
  name="${entry%%|*}"; path="${entry#*|}"
  if [ ! -f "$path" ]; then echo "SKIP (missing): $name  $path"; continue; fi
  echo "================ $name ================"
  python 00_compute_metrics.py --latent "$path" --method "$name" --label "$LABEL" --out .
done

python 01_adjust_accuracy.py --label "$LABEL" --out .

conda deactivate
echo "DONE: sum_metrics.csv + celltype_metrics.csv + adj_atac_predaccu.csv (SingleR-labelled) -> run rank_stability.py"
