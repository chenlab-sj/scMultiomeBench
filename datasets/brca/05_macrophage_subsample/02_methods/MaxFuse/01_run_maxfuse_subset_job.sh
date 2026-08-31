#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

# One MaxFuse (MACROSUB, SUB) job: Signac prep (seurat4 env) -> Fusor (maxfuse-env). Reads MACROSUB + SUB
# from env. Both stages write only under {MACROSUB}/maxfuse/sub{SUB}/ -> isolated, parallel-safe.
set -u
SUBS=${PROJECT_ROOT}/BRCA/HT243-S1H4_subsample
O="$SUBS/$MACROSUB/maxfuse/sub$SUB"
mkdir -p "$O/input"

echo "=== [1/2] Signac prep: $MACROSUB sub$SUB (seurat4) ==="
module load conda3/202210
conda activate seurat4
R CMD BATCH --no-save "$SUBS/prep_maxfuse_subset.R" "$O/prep.Rout"
conda deactivate
module unload conda3/202210 2>/dev/null

# fail early if the prep didn't produce all 4 inputs (R CMD BATCH returns 0 even on R error)
for h in maxfuse_rna_shared maxfuse_atac_shared maxfuse_rna_active maxfuse_atac_active; do
  [ -f "$O/input/$h.h5ad" ] || { echo "PREP FAILED: missing $h.h5ad (see $O/prep.Rout)"; exit 1; }
done

echo "=== [2/2] Fusor: $MACROSUB sub$SUB (maxfuse-env) ==="
module load conda3/202311
source activate maxfuse-env
export PYTHONNOUSERSITE=1
python "$SUBS/run_maxfuse_subset.py"
conda deactivate
echo "DONE: $O/latent.csv"
