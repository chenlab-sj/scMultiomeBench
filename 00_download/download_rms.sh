#!/usr/bin/env bash
#
# scMultiomeBench 00_download :: RMS Mast607A (+ two training samples, open access)
#
# The RMS arm starts from the authors' own 10x Multiome sequencing, processed
# with cellranger-arc count v2.0.0 against GRCh37/hg19 (see
# datasets/rms/01_preprocess/00_cellranger_arc_sub.sh). This script documents
# the route and verifies the expected layout; the cellranger-arc step is
# site-specific and is not run from here.
#
# GEO accession: GSE209784 (open access).
#
# Expected layout, as read by datasets/rms/01_preprocess/ (note: RMS/, not rms/):
#   $DATA_ROOT/RMS/<sample>/<sample>/outs/filtered_feature_bc_matrix.h5
#   $DATA_ROOT/RMS/<sample>/<sample>/outs/atac_fragments.tsv.gz (+ .tbi)
#   $DATA_ROOT/RMS/<sample>/<sample>/outs/filtered_feature_bc_matrix/features.tsv
# for samples Mast607A_TB19_22652 (test), Mast39_TB12_1442 and
# MAST213F_TB15_5705 (training/reference).
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../config/config.sh
. "${HERE}/../config/config.sh"

case "${DATA_ROOT}" in /path/to/data*)
  echo "REFUSING to run: DATA_ROOT is still the placeholder ${DATA_ROOT}." >&2
  echo "Create config/config.local.sh from config/config.local.sh.example first." >&2
  exit 2 ;;
esac

BASE="${DATA_ROOT}/RMS"
missing=0
for sample in Mast607A_TB19_22652 Mast39_TB12_1442 MAST213F_TB15_5705; do
  for f in "outs/filtered_feature_bc_matrix.h5" "outs/atac_fragments.tsv.gz" "outs/atac_fragments.tsv.gz.tbi"; do
    p="${BASE}/${sample}/${sample}/${f}"
    if [ -s "${p}" ]; then echo "have    ${sample}/${f}"
    else echo "MISSING ${sample}/${f}"; missing=1; fi
  done
done

[ "${missing}" = 0 ] && { echo "done: RMS layout complete under ${BASE}"; exit 0; }

cat >&2 <<EOF

MANUAL STEPS REQUIRED (open access, but compute-heavy)
  1. The GEO accession is GSE209784.
  2. Download the raw sequencing for the three samples from GEO/SRA.
  3. Process each sample with cellranger-arc count v2.0.0 against GRCh37/hg19.
  4. Place each sample's Cell Ranger 'outs/' as shown above under
     ${BASE}/<sample>/<sample>/
  Then run datasets/rms/01_preprocess/01_make_common_peak_h5.R to build the
  common-peak matrices.
EOF
exit 1
