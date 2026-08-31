#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J figS1c_metrics
#BSUB -q large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o figS1c_metrics.log
#BSUB -e figS1c_metrics.err

# FigS1C (pbmc10k composite table) -- FULL recompute of all 23 methods with ONE pipeline
# (00_compute_metrics.py, same scib/KS defs as the original benchmark_metrics.py), then the
# random-adjusted ATAC accuracy. 18 original latents live on ${TOOLS_ROOT}/...; the 5 new
# ones under scripts/<m>/latent.csv. Runs in benchmark_env (scib + benchmark_fun). Peak score
# is added separately (02_peak_similarity_sub.sh: old copied peak + new-method peak).
module load conda3/202303
conda activate benchmark_env

ROOT=${PROJECT_ROOT}
cd "${ROOT}/pbmc/pbmc10k/benchmark"
export BENCHMARK_FUN_DIR=${BENCHMARK_FUN_DIR}
LABEL=label.csv

# name|latent_path  (18 original = the paths benchmark_metrics.py used; 5 new = repo)
methods=(
  "BindSC|${TOOLS_ROOT}/bindsc/pbmc10k/res_pbmc10k_control/bindsc_lat_df.csv"
  "Conos|${TOOLS_ROOT}/Concos/pbmc10k/res_pbmc10k_control/coembed_coor.csv"
  "LIGER|${TOOLS_ROOT}/liger/pbmc10k/res_pbmc10k_control/coembed_coor.csv"
  "MultiMAP|${TOOLS_ROOT}/multimap/pbmc10k/res_pbmc10k_control/coembed_coor.csv"
  "Seurat(CCA)|${TOOLS_ROOT}/Seuratv3/pbmc10k/res_pbmc10k_control/lat_df.csv"
  "scDART|${TOOLS_ROOT}/scDART/scDART/pbmc10k/res_pbmc10k/latent.csv"
  "scMoMaT|${TOOLS_ROOT}/scMoMaT/pbmc10k/latent.csv"
  "simba|${TOOLS_ROOT}/simba/pbmc10k/latent.csv"
  "scglue|${TOOLS_ROOT}/scglue/pbmc10k/res_pbmc10K_control/latent.csv"
  "Unioncom|${TOOLS_ROOT}/unioncom/pbmc10k/latent.csv"
  "Portal|${TOOLS_ROOT}/portal/Portal/pbmc10k/lat_df.csv"
  "scJoint|${TOOLS_ROOT}/scJoint/pbmc10k/res_control/scjoint_lat_df.csv"
  "scBridge|${TOOLS_ROOT}/scBridge/pbmc10k/latent.csv"
  "Cobolt|${TOOLS_ROOT}/cobolt/pbmc10k/res_pbmc10k_test_all/cobolt_latent.csv"
  "scglue(multiome)|${TOOLS_ROOT}/scglue/pbmc10k/res_pbmc10K_withpair_control/scglue_latent.csv"
  "scVI|${TOOLS_ROOT}/scvi/pbmc10k/res_pbmc10k_test_all/scvi_latent.csv"
  "Seurat(WNN)|${TOOLS_ROOT}/Seuratv4/pbmc10k/res_pbmc10k_test_all/seurat4_latent.csv"
  "MinNet|${TOOLS_ROOT}/MinNet/pbmc10k_control/MinNet_lat_df.csv"
  # ---- new methods (this revision) ----
  "MaxFuse|${ROOT}/pbmc/pbmc10k/scripts/maxfuse/latent.csv"
  "Multigrate|${ROOT}/pbmc/pbmc10k/scripts/multigrate/latent.csv"
  "MIDAS|${ROOT}/pbmc/pbmc10k/scripts/midas/latent.csv"
  "scButterfly|${ROOT}/pbmc/pbmc10k/scripts/scbutterfly/latent.csv"
  "MIRA|${ROOT}/pbmc/pbmc10k/scripts/MIRA/latent.csv"
)

for entry in "${methods[@]}"; do
  name="${entry%%|*}"; path="${entry#*|}"
  if [ ! -f "$path" ]; then echo "SKIP (missing): $name  $path"; continue; fi
  echo "================ $name ================"
  python 00_compute_metrics.py --latent "$path" --method "$name" --label "$LABEL" --out .
done

# random-adjusted ATAC accuracy (average_accu) from all per-method predicted labels
python 01_adjust_accuracy.py --label "$LABEL" --out .

conda deactivate
echo "DONE: sum_metrics.csv, celltype_metrics.csv, knn_pred_accu.csv, adj_atac_predaccu.csv"
