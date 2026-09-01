#!/usr/bin/env bash
#
# multiomeBench 00_download :: BMMC s1d1 / s2d1 / s4d1 (open access)
#
# The BMMC arm cannot be fetched as ready-made files: the authors obtained raw
# sequencing from NCBI SRA and reprocessed it with Cell Ranger ARC, and the
# preprocessing scripts read the resulting Cell Ranger output layout. This
# script therefore documents the route and verifies the expected layout; it
# does not download.
#
# GEO accession: GSE194122 (NeurIPS 2021 Open Problems 10x Multiome BMMC;
# the s1d1/s2d1/s4d1 site-donor naming follows that benchmark). Confirmed by
# the authors -- see DATA.md section 6.
#
# Expected layout, as read by datasets/bmmc_d1/01_preprocess/00_file_prep.R
# (note: BMMC/NCBI_sra/, not bmmc_d1/):
#   $DATA_ROOT/BMMC/NCBI_sra/<batch>/outs/filtered_feature_bc_matrix.h5
#   $DATA_ROOT/BMMC/NCBI_sra/<batch>/outs/atac_fragments.tsv.gz (+ .tbi)
# for batches s1d1, s2d1, s4d1.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../config/config.sh
. "${HERE}/../config/config.sh"

case "${DATA_ROOT}" in /path/to/data*)
  echo "REFUSING to run: DATA_ROOT is still the placeholder ${DATA_ROOT}." >&2
  echo "Create config/config.local.sh from config/config.local.sh.example first." >&2
  exit 2 ;;
esac

BASE="${DATA_ROOT}/BMMC/NCBI_sra"
missing=0
for batch in s1d1 s2d1 s4d1; do
  for f in "outs/filtered_feature_bc_matrix.h5" "outs/atac_fragments.tsv.gz" "outs/atac_fragments.tsv.gz.tbi"; do
    p="${BASE}/${batch}/${f}"
    if [ -s "${p}" ]; then echo "have    ${batch}/${f}"
    else echo "MISSING ${batch}/${f}"; missing=1; fi
  done
done

[ "${missing}" = 0 ] && { echo "done: BMMC layout complete under ${BASE}"; exit 0; }

cat >&2 <<EOF

MANUAL STEPS REQUIRED (open access, but compute-heavy)
  1. The GEO accession is GSE194122 (NeurIPS 2021 Open Problems 10x Multiome
     BMMC); identify the SRA runs for site-donor batches s1d1, s2d1, s4d1.
  2. Download the raw sequencing for those runs from SRA
     (prefetch/fasterq-dump, or the SRA cloud mirrors).
  3. Process each batch with cellranger-arc count against GRCh38.
  4. Place each batch's Cell Ranger 'outs/' as shown above under
     ${BASE}/<batch>/
EOF
exit 1
