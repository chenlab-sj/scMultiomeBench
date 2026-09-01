#!/usr/bin/env bash
#
# scMultiomeBench 00_download :: BRCA HT243B1-S1H4 (CONTROLLED ACCESS -- dbGaP)
#
# This script CANNOT download anything: the BRCA data requires an approved
# dbGaP Data Access Request (study phs002371.v3.p1, HTAN WUSTL) followed by an
# authenticated Synapse download, neither of which can be automated in a
# public script. It verifies the expected layout and prints the access route.
# See DATA.md section 7 for the full instructions.
#
# Expected layout, as read by datasets/brca/01_preprocess/ (note: HTAN/, the
# original working-directory name, not brca/):
#   $DATA_ROOT/HTAN/HT243B1-S1H4/syn53214720/HT243B1-S1H4.rds   (RNA, test)
#   $DATA_ROOT/HTAN/HT243B1-S1H4/syn53215789/HT243B1-S1H4.rds   (ATAC, test)
#   $DATA_ROOT/HTAN/HT243B1-S1H4/HT243B1-S1H4-atac_fragments.tsv.gz (+ .tbi, tabix -p bed)
#   $DATA_ROOT/HTAN/HT263B1-S1H1/syn53214683/HT263B1-S1H1.rds   (RNA, paired training)
#   $DATA_ROOT/HTAN/HT263B1-S1H1/syn53215774/HT263B1-S1H1.rds   (ATAC, paired training)
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../config/config.sh
. "${HERE}/../config/config.sh"

case "${DATA_ROOT}" in /path/to/data*)
  echo "REFUSING to run: DATA_ROOT is still the placeholder ${DATA_ROOT}." >&2
  echo "Create config/config.local.sh from config/config.local.sh.example first." >&2
  exit 2 ;;
esac

BASE="${DATA_ROOT}/HTAN"
missing=0
for f in \
  "HT243B1-S1H4/syn53214720/HT243B1-S1H4.rds" \
  "HT243B1-S1H4/syn53215789/HT243B1-S1H4.rds" \
  "HT243B1-S1H4/HT243B1-S1H4-atac_fragments.tsv.gz" \
  "HT243B1-S1H4/HT243B1-S1H4-atac_fragments.tsv.gz.tbi" \
  "HT263B1-S1H1/syn53214683/HT263B1-S1H1.rds" \
  "HT263B1-S1H1/syn53215774/HT263B1-S1H1.rds"; do
  p="${BASE}/${f}"
  if [ -s "${p}" ]; then echo "have    ${f}"
  else echo "MISSING ${f}"; missing=1; fi
done

[ "${missing}" = 0 ] && { echo "done: BRCA layout complete under ${BASE}"; exit 0; }

cat >&2 <<EOF

MANUAL STEPS REQUIRED (CONTROLLED ACCESS -- cannot be scripted)
  1. Obtain dbGaP authorization for study phs002371.v3.p1 (HTAN WUSTL):
     submit a Data Access Request at https://dbgap.ncbi.nlm.nih.gov/
     (institutional signing official, research use statement; turnaround is
     typically weeks).
  2. Once approved, download the four Seurat objects from Synapse
     (https://www.synapse.org/, account linked to the approved credentials):
       synapse get syn53214720 syn53215789   # HT243B1-S1H4 RNA / ATAC
       synapse get syn53214683 syn53215774   # HT263B1-S1H1 RNA / ATAC (training)
  3. Obtain the HT243B1-S1H4 ATAC fragments file and index it:
     tabix -p bed HT243B1-S1H4-atac_fragments.tsv.gz
  4. Place everything as shown in this script's header under ${BASE}/
  Do NOT redistribute any BRCA object, fragment file, barcode list, or
  per-cell table. See DATA.md section 7.
EOF
exit 1
