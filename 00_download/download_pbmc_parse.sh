#!/usr/bin/env bash
#
# multiomeBench 00_download :: Parse PBMC (Parse Biosciences, open access)
#
# The Parse Evercode WT Mini v3 PBMC dataset is distributed from a download
# page rather than a stable direct URL, so this script cannot hard-code a
# link that is guaranteed to stay valid. It verifies the expected layout and,
# if PARSE_DGE_URL is set to the direct link of the DGE archive from the page
# below, fetches and unpacks it.
#
#   https://www.parsebiosciences.com/datasets/performance-of-evercode-wt-mini-v3-in-human-pbmcs/#download
#
# Files consumed by datasets/pbmc_parse/01_preprocess/00_data_prep.R
# (via Seurat::ReadParseBio): DGE.mtx, all_genes.csv, cell_metadata.csv.
# The ATAC side of this benchmark is reused from pbmc3k and needs no download.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../config/config.sh
. "${HERE}/../config/config.sh"

case "${DATA_ROOT}" in /path/to/data*)
  echo "REFUSING to run: DATA_ROOT is still the placeholder ${DATA_ROOT}." >&2
  echo "Create config/config.local.sh from config/config.local.sh.example first." >&2
  exit 2 ;;
esac

DEST="${DATA_ROOT}/pbmc_parse"
mkdir -p "${DEST}"

have_all() {
  [ -s "${DEST}/DGE.mtx" ] && [ -s "${DEST}/all_genes.csv" ] && [ -s "${DEST}/cell_metadata.csv" ]
}

if have_all; then
  echo "have    Parse DGE triplet in ${DEST} -- nothing to do"
  exit 0
fi

if [ -n "${PARSE_DGE_URL:-}" ]; then
  echo "fetch   ${PARSE_DGE_URL}"
  archive="${DEST}/parse_dge_download"
  curl -fL --retry 3 -C - -o "${archive}" "${PARSE_DGE_URL}"
  case "$(file -b "${archive}" 2>/dev/null)" in
    Zip*)  unzip -o -j "${archive}" -d "${DEST}" ;;
    gzip*) tar -xzf "${archive}" -C "${DEST}" --strip-components 1 2>/dev/null \
             || tar -xzf "${archive}" -C "${DEST}" ;;
    *)     echo "unrecognised archive type; unpack ${archive} by hand into ${DEST}" >&2 ;;
  esac
  have_all && { echo "done: Parse DGE triplet ready in ${DEST}"; exit 0; }
fi

cat >&2 <<EOF
MANUAL STEP REQUIRED
  1. Open the download page (open access, no registration):
     https://www.parsebiosciences.com/datasets/performance-of-evercode-wt-mini-v3-in-human-pbmcs/#download
  2. Download the DGE (digital gene expression) output.
  3. Place DGE.mtx, all_genes.csv and cell_metadata.csv in:
     ${DEST}
  (or rerun this script with PARSE_DGE_URL=<direct link> to fetch it here)
EOF
exit 1
