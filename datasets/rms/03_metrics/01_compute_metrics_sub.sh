#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J rms_figS4a_metrics
#BSUB -q large_mem
#BSUB -R rusage[mem=120000]
#BSUB -n 1
#BSUB -o figS4a_metrics.log
#BSUB -e figS4a_metrics.err

# RMS Mast607A benchmark matrix (FigS4A) -- same parameterized 00_compute_metrics.py + 01_adjust_accuracy.py
# as pbmc3k/BRCA, on the RMS label.csv. One compute_metrics call per method on its rep1 (base) latent,
# read from the staged dir (../staged/<method>/rep1.csv, normalized by 01_prep_latents.py).
# Missing/format-bad latents are SKIPPED (loop continues) -> simba/scButterfly fill in once they finish.
module load conda3/202303
conda activate benchmark_env

ROOT=${PROJECT_ROOT}
cd "${ROOT}/RMS/benchmark/figS4a"
export BENCHMARK_FUN_DIR=${BENCHMARK_FUN_DIR}
LABEL=label.csv
STAGE=../staged

python ../prep_latents.py        # (re)stage all 14 methods x 3 reps

methods=(scVI Cobolt Portal scglue "scglue(multiome)" BindSC "Seurat(CCA)" \
         MaxFuse MIDAS scDART scBridge simba scButterfly scJoint)
for name in "${methods[@]}"; do
  path="${STAGE}/${name}/rep1.csv"
  if [ ! -f "$path" ]; then echo "SKIP (missing): $name  $path"; continue; fi
  echo "================ $name ================"
  python 00_compute_metrics.py --latent "$path" --method "$name" --label "$LABEL" --out . || echo "FAILED: $name"
done

python 01_adjust_accuracy.py --label "$LABEL" --out .

conda deactivate
echo "DONE: sum_metrics.csv, celltype_metrics.csv, knn_pred_accu.csv, adj_atac_predaccu.csv"
echo "      (plot separately: Rscript plot_metrics_matrix.R)"
