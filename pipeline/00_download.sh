#!/usr/bin/env bash
#
# multiomeBench -- pipeline stage 0: obtain the raw datasets
#
# WHAT THIS IS
#   A thin driver that records where each benchmark dataset comes from and, if
#   fetch scripts are present in 00_download/, the order in which to run them.
#   It contains no analysis logic of its own.
#
# DRY RUN IS THE DEFAULT.
#   With DRY_RUN=1 (the default) this script only PRINTS what it would do.
#   Nothing is downloaded and nothing is written.
#   To actually execute:   DRY_RUN=0 ./00_download.sh [dataset]
#
# USAGE
#   ./00_download.sh                    # dry run, all datasets
#   ./00_download.sh rms                # dry run, one dataset
#   DRY_RUN=0 ./00_download.sh pbmc3k   # actually execute, one dataset
#
#   datasets: pbmc3k pbmc10k pbmc_parse bmmc_d1 brca rms
#
# STATE OF THIS STAGE
#   00_download/ holds one fetch script per open dataset (download_<dataset>.sh).
#   The 10x PBMC scripts download directly; the Parse, BMMC and RMS scripts
#   verify the expected layout and document the manual route where a stable
#   direct URL or a scriptable pipeline does not exist. This driver runs them
#   in filename order and then prints the accession notes below.
#
# PATHS
#   Sourced from ../config/config.sh, which ships PLACEHOLDER defaults
#   (DATA_ROOT=/path/to/data). Copy config/config.local.sh.example to
#   config/config.local.sh and set real paths before DRY_RUN=0.
#
# MARKERS in the printed output
#   [dry-run]  a command that would be executed
#   [manual]   a step that has no script; do it by hand from the stated accession
#   #          a note
#
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../config/config.sh
. "${HERE}/../config/config.sh"

: "${DRY_RUN:=1}"
WANT="${1:-all}"

if [ "${DRY_RUN}" != "1" ]; then
  case "${DATA_ROOT}" in
    /path/to/data*)
      echo "REFUSING to run: DATA_ROOT is still the placeholder ${DATA_ROOT}." >&2
      echo "Create config/config.local.sh from config/config.local.sh.example first." >&2
      exit 2 ;;
  esac
  echo "### DRY_RUN=0 -- commands WILL be executed. DATA_ROOT=${DATA_ROOT}"
fi

want()   { case "${WANT}" in all|"$1") return 0 ;; *) return 1 ;; esac }
say()    { printf '\n=== %s\n' "$*" ; }
note()   { printf '  # %s\n' "$*" ; }
manual() { printf '  [manual]   %s\n' "$*" ; }

run() {
  local rel="$1"; shift
  local abs="${PROJECT_ROOT}/${rel}"
  local -a cmd
  case "${rel}" in
    *.R)  cmd=(Rscript "$(basename "${rel}")") ;;
    *.py) cmd=(python  "$(basename "${rel}")") ;;
    *.sh) cmd=(bash    "$(basename "${rel}")") ;;
    *)    printf '  [skip]     no interpreter for %s\n' "${rel}"; return 0 ;;
  esac
  if [ "${DRY_RUN}" = "1" ]; then
    printf '  [dry-run]  (cd %s && %s)\n' "$(dirname "${rel}")" "${cmd[*]}${*:+ $*}"
  else
    if [ ! -f "${abs}" ]; then printf '  MISSING    %s\n' "${rel}" >&2; return 1; fi
    ( cd "$(dirname "${abs}")" && "${cmd[@]}" "$@" )
  fi
}

# Any scripts the authors (or you) add to 00_download/ are run first, in filename
# order, so the numeric prefixes 00_/01_/... control the sequence.
say "00 download :: scripts in 00_download/"
_found=0
for f in "${PROJECT_ROOT}"/00_download/*; do
  [ -e "$f" ] || continue
  _found=1
  run "00_download/$(basename "$f")"
done
[ "${_found}" = 0 ] && note "00_download/ contains no scripts; every dataset below is fetched by hand"

# --------------------------------------------------------------------------
# Accessions. Everything except BRCA is openly available.
# --------------------------------------------------------------------------

if want pbmc3k; then
  say "00 download :: PBMC 3k"
  note "10x Genomics public single-cell multiome (RNA + ATAC) demonstration dataset"
  note "open access; no application required"
  manual "download the filtered feature-barcode matrix, the ATAC fragments file and its index"
  manual "place them under \${DATA_ROOT}/pbmc3k/  (DATA_ROOT=${DATA_ROOT})"
fi

if want pbmc10k; then
  say "00 download :: PBMC 10k"
  note "10x Genomics public single-cell multiome (RNA + ATAC) demonstration dataset"
  note "open access; no application required"
  manual "download the filtered feature-barcode matrix, the ATAC fragments file and its index"
  manual "place them under \${DATA_ROOT}/pbmc10k/"
fi

if want pbmc_parse; then
  say "00 download :: Parse PBMC"
  note "Parse Biosciences Evercode WT Mini v3 public PBMC dataset"
  note "open access; no application required"
  note "used as the held-out generalisation dataset for Fig6 / FigS12"
  manual "download the Parse count matrices and place them under \${DATA_ROOT}/pbmc_parse/"
fi

if want bmmc_d1; then
  say "00 download :: BMMC s1d1"
  note "public bone-marrow mononuclear cell single-cell multiome dataset"
  note "GEO accession GSE194122; open access, no application required"
  note "site s1 donor d1 is the primary sample; other sites are used for the cross-site analysis"
  manual "download the multiome matrices and place them under \${DATA_ROOT}/bmmc_d1/"
fi

if want brca; then
  say "00 download :: BRCA (CONTROLLED ACCESS)"
  note "HTAN WUSTL breast cancer; dbGaP accession phs002371.v3.p1"
  note "THIS IS THE ONLY DATASET THAT IS NOT OPEN. It requires an approved dbGaP"
  note "data access request. No script here can fetch it for you."
  note "Samples used:"
  note "  HT243B1-S1H4  -- the benchmarked sample.  Synapse syn53214720 (RNA) / syn53215789 (ATAC)"
  note "  HT263B1-S1H1  -- the paired multiome TRAINING sample for HT243."
  note "                   Synapse syn53214683 (RNA) / syn53215774 (ATAC)"
  manual "obtain dbGaP approval for phs002371.v3.p1, then download the four Synapse objects above"
  manual "place them under \${DATA_ROOT}/brca/"
fi

if want rms; then
  say "00 download :: RMS Mast607A"
  note "rhabdomyosarcoma single-cell multiome generated by the authors of this benchmark"
  note "GEO accession GSE209784; open access"
  manual "download GSE209784 and place the combined-multiome h5 and fragments under \${DATA_ROOT}/rms/"
  note "stage 01 re-derives a common peak set for this dataset (01_make_common_peak_h5.R);"
  note "the cellranger-arc step that precedes it is datasets/rms/01_preprocess/00_cellranger_arc_sub.sh,"
  note "an LSF submitter that is site-specific and is not run by this driver"
fi

say "stage 0 listing complete (DRY_RUN=${DRY_RUN})"
