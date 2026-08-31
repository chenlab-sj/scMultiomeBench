#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

# Copy the scglue / scglue_withpair / portal rep latents from their cluster run folders into the repo,
# so all three reps (rep1=base latent, rep2, rep3) live under RMS/Mast607/script/<method>/{,rep2,rep3}/.
# These 3 methods write to ${TOOLS_ROOT}/...; everything else already writes into the repo.
# Run this ON THE CLUSTER (it reads ${TOOLS_ROOT}/... and writes the mounted /path/to/project/... repo).

REPO=${PROJECT_ROOT}/RMS/Mast607/script

copy_one () {
  local m=$1 base=$2 fn=$3
  echo "==== $m ===="
  mkdir -p "$REPO/$m"
  # rep1 = the base latent (the existing benchmark run)
  if cp "$base/$fn" "$REPO/$m/$fn" 2>/dev/null; then echo "  rep1  <- $base/$fn"; else echo "  rep1  MISSING: $base/$fn"; fi
  # rep2 / rep3
  for r in rep2 rep3; do
    mkdir -p "$REPO/$m/$r"
    if cp "$base/$r/$fn" "$REPO/$m/$r/$fn" 2>/dev/null; then echo "  $r  <- $base/$r/$fn"; else echo "  $r  MISSING: $base/$r/$fn"; fi
  done
}

copy_one scglue          ${TOOLS_ROOT}/scglue/RMS/Mast607A_TB19_22652                  latent.csv
copy_one scglue_withpair ${TOOLS_ROOT}/scglue/scglue_withpair/RMS/Mast607A_TB19_22652  scglue_latent.csv
copy_one portal          ${DATA_ROOT}/RMS/Mast607A_TB19_22652/portal             lat_df.csv

echo ""
echo "DONE. Any 'MISSING' line above = that rep hasn't finished yet (re-run it, then re-run this copy)."
