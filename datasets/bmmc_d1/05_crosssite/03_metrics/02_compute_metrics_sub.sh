#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J bmmc_xs_metrics
#BSUB -q large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o cs_metrics.log
#BSUB -e cs_metrics.err

# BMMC same-donor cross-site benchmark matrix (step 1). Builds label.csv + label_major.csv, normalizes
# every method latent to bare <bc>-1_rna / <bc>-1_atac (01_prep_latents.py), then runs 00_compute_metrics.py
# per method on label.csv -> sum/celltype/accu + knn_pred_label__*.csv (the peak step reads those).
# GENUINELY UNPAIRED (crosssite) / paired (s1d1_paired): the same-cell metrics are dropped at plot time.
# LS_SUBCWD = the experiment dir this was `bsub <`-submitted from, so ONE script serves both experiments.
module load conda3/202303
conda activate benchmark_env

FIG="${LS_SUBCWD:-$PWD}"          # .../BMMC_d1/benchmark/{crosssite,s1d1_paired}
cd "${FIG}" || exit 1
export BENCHMARK_FUN_DIR=${BENCHMARK_FUN_DIR}
LABEL=label.csv

# ---- build label.csv + label_major.csv, normalize all latents -> latents/<key>.csv ----
python 00_build_label.py
python 01_prep_latents.py

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
  if [ ! -f "$path" ]; then echo "SKIP (missing): $name  $path"; continue; fi
  echo "================ $name ================"
  python 00_compute_metrics.py --latent "$path" --method "$name" --label "$LABEL" --out . || echo "FAILED: $name"
done

conda deactivate
echo "DONE: sum_metrics.csv, celltype_metrics.csv, knn_pred_label__*.csv (then: major + peak -> plot)"
