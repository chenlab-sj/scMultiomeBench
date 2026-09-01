#!/usr/bin/env bash
#
# multiomeBench 00_download :: PBMC 3k (10x Genomics, open access)
#
# Fetches the inputs that datasets/pbmc3k/01_preprocess/ reads, into
# $DATA_ROOT/pbmc3k/, keeping the provider's filenames. Idempotent: a file
# already present (and passing gzip integrity where applicable) is not refetched.
#
# Source: 10x Genomics dataset portal, "PBMC from a healthy donor -
# granulocytes removed through cell sorting (3k)", Cell Ranger ARC 1.0.0.
# Verify the md5 sums shown on the dataset page after download.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../config/config.sh
. "${HERE}/../config/config.sh"

case "${DATA_ROOT}" in /path/to/data*)
  echo "REFUSING to run: DATA_ROOT is still the placeholder ${DATA_ROOT}." >&2
  echo "Create config/config.local.sh from config/config.local.sh.example first." >&2
  exit 2 ;;
esac

BASE="https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_3k"
DEST="${DATA_ROOT}/pbmc3k"
mkdir -p "${DEST}"

fetch() {
  local name="$1" out="${DEST}/$1"
  if [ -s "${out}" ]; then
    case "${name}" in *.gz) gzip -t "${out}" 2>/dev/null && { echo "have    ${name}"; return 0; } ;;
                      *)    echo "have    ${name}"; return 0 ;; esac
    echo "corrupt ${name} -- refetching"
  fi
  echo "fetch   ${name}"
  curl -fL --retry 3 -C - -o "${out}.part" "${BASE}/${name}" && mv "${out}.part" "${out}"
}

fetch pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5
fetch pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz
fetch pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz.tbi

echo "done: $(ls "${DEST}" | wc -l | tr -d ' ') files in ${DEST}"
