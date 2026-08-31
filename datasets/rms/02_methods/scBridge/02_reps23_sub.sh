#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"


#BSUB -P benchmark
#BSUB -J rms_scbridge_reps
#BSUB -q dgx
#BSUB -gpu "num=1"
#BSUB -R rusage[mem=150000]
#BSUB -o run_scbridge_reps.log
#BSUB -e run_scbridge_reps.err

# scBridge rep2 + rep3 (seeds 7, 42). scBridge is a TWO-STEP, CLONED-repo pipeline run from
# ${TOOLS_ROOT}/scBridge:
#   s1) python main.py --data_path=<DATA_DIR> --source_data=scBridge_rna.h5ad --target_data=scBridge_atac.h5ad --random_seed=<SEED> --umap_plot
#       -> trains, dumps tmp_saved_objects.pkl (in repo dir) + *-integrated.h5ad into <DATA_DIR>
#   s2) python 2_save_scBridge_result.py
#       -> reads tmp_saved_objects.pkl, writes latent.csv + label_pred.csv into args.data_path (= <DATA_DIR>)
#
# The ONLY per-rep knob is --random_seed (this is exactly how BRCA + pbmc reps were produced:
# rep dirs differ from base only by --random_seed; data_path stays the same shared dir).
# 2_save_scBridge_result.py has NO output-path arg -- it always overwrites latent.csv/label_pred.csv
# in <DATA_DIR>. So the reps MUST be serialized: per seed, run s1 then s2, then COPY the freshly
# written latent.csv/label_pred.csv into the rep dir before the next seed overwrites them.
#
# rep1 = the existing latent.csv in DATA_DIR (untouched). 7/42 avoid the prior rep seeds (0/40)
# and the base run (default seed). main.py uses the GPU -> dgx.

module load conda3/202210
conda activate scBridge
module load gcc/13.1.0-rhel7

REPO=${TOOLS_ROOT}/scBridge
# repo dir: the scBridge_rna/atac.h5ad inputs + the base latent.csv live HERE (not ${TOOLS_ROOT}/scBridge/RMS/)
DATA_DIR=${PROJECT_ROOT}/RMS/Mast607/script/scBridge

cd "$REPO"

# preserve the existing base latent (rep1) before the reps overwrite DATA_DIR/latent.csv
mkdir -p "$DATA_DIR/rep1"
cp "$DATA_DIR/latent.csv"     "$DATA_DIR/rep1/latent.csv"     2>/dev/null || true
cp "$DATA_DIR/label_pred.csv" "$DATA_DIR/rep1/label_pred.csv" 2>/dev/null || true

for rs in "rep2 7" "rep3 42"; do
  set -- $rs; REP=$1; SEED=$2
  echo "=== scBridge $REP (seed $SEED) ==="

  # main.py builds the input path by STRING-CONCAT (data_path + source_data) with NO separator, so
  # DATA_DIR MUST end in '/' -- otherwise it looks for .../scBridgescBridge_rna.h5ad and FileNotFounds,
  # main.py dies, and 2_save silently re-emits the STALE tmp_saved_objects.pkl (-> identical reps).
  # 2_save reads the *-integrated.h5ad that main.py writes -- clear the prior rep's so a failed train
  # can't make 2_save silently re-emit the previous rep's latent.
  rm -f "$DATA_DIR"/*-integrated.h5ad "$REPO"/tmp_saved_objects.pkl

  # s1: train (writes *-integrated.h5ad into DATA_DIR + tmp_saved_objects.pkl in REPO)
  python main.py --data_path="$DATA_DIR/" --source_data="scBridge_rna.h5ad" \
      --target_data="scBridge_atac.h5ad" --random_seed="$SEED" --umap_plot

  # s2: 2_save ALSO concats args.data_path + filenames and needs source/target -- pass the SAME args
  # (called bare it errors on None+str and the cp below grabs the stale latent -> identical reps).
  python 2_save_scBridge_result.py --data_path="$DATA_DIR/" \
      --source_data="scBridge_rna.h5ad" --target_data="scBridge_atac.h5ad"

  # stash this rep's outputs before the next seed overwrites the shared files
  mkdir -p "$DATA_DIR/$REP"
  cp "$DATA_DIR/latent.csv"     "$DATA_DIR/$REP/latent.csv"
  cp "$DATA_DIR/label_pred.csv" "$DATA_DIR/$REP/label_pred.csv"
  echo "=== $REP done -> $DATA_DIR/$REP/latent.csv ==="
done

conda deactivate
echo "DONE: $DATA_DIR/rep2,rep3/latent.csv"
