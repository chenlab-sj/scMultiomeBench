#!/usr/bin/env bash
#
# multiomeBench -- pipeline stage 1: per-dataset preprocessing
#
# WHAT THIS IS
#   A thin driver that records the execution ORDER of the real analysis scripts.
#   It contains no analysis logic of its own; every line below points at a script
#   that lives under datasets/ (or summary/) and does the actual work.
#
# DRY RUN IS THE DEFAULT.
#   With DRY_RUN=1 (the default) this script only PRINTS the commands it would run.
#   Nothing is executed, no file is written, no job is submitted.
#   To actually execute:   DRY_RUN=0 ./01_preprocess.sh [dataset]
#   Do not do that casually -- see pipeline/README.md for the runtime cost.
#
# USAGE
#   ./01_preprocess.sh                 # dry run, all datasets
#   ./01_preprocess.sh pbmc3k          # dry run, one dataset
#   DRY_RUN=0 ./01_preprocess.sh rms   # actually execute, one dataset
#
#   datasets: pbmc3k pbmc10k pbmc_parse bmmc_d1 brca rms   (and "summary" where present)
#
# PATHS
#   Sourced from ../config/config.sh, which ships PLACEHOLDER defaults
#   (DATA_ROOT=/path/to/data). Copy config/config.local.sh.example to
#   config/config.local.sh and set real paths before DRY_RUN=0.
#
# MARKERS in the printed output
#   [dry-run]  a command that would be executed
#   [lsf]      an original LSF submission script (bsub); site-specific, never auto-run
#   [config]   a file you edit rather than execute (e.g. scJoint 02_config.py)
#   [lib]      an imported library, not an entry point
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

want() { case "${WANT}" in all|"$1") return 0 ;; *) return 1 ;; esac }
say()  { printf '\n=== %s\n' "$*" ; }
sub()  { printf -- '\n--- %s\n' "$*" ; }
note() { printf '  # %s\n' "$*" ; }
lsf()  { printf '  [lsf]      %s\n' "$1" ; }
cfg()  { printf '  [config]   %s\n' "$1" ; }
lib()  { printf '  [lib]      %s\n' "$1" ; }

run() {
  local rel="$1"; shift
  local abs="${PROJECT_ROOT}/${rel}"
  local -a cmd
  case "${rel}" in
    *.R)     cmd=(Rscript "$(basename "${rel}")") ;;
    *.py)    cmd=(python  "$(basename "${rel}")") ;;
    *.ipynb) cmd=(jupyter nbconvert --to notebook --execute --inplace "$(basename "${rel}")") ;;
    *.sh)    cmd=(bash    "$(basename "${rel}")") ;;
    *)       printf '  [skip]     no interpreter for %s\n' "${rel}"; return 0 ;;
  esac
  if [ "${DRY_RUN}" = "1" ]; then
    printf '  [dry-run]  (cd %s && %s)\n' "$(dirname "${rel}")" "${cmd[*]}${*:+ $*}"
  else
    if [ ! -f "${abs}" ]; then printf '  MISSING    %s\n' "${rel}" >&2; return 1; fi
    ( cd "$(dirname "${abs}")" && "${cmd[@]}" "$@" )
  fi
}

# Consumes: raw matrices / fragments under $DATA_ROOT (see 00_download.sh).
# Produces: the per-dataset RNA and ATAC objects, the common-peak matrices and the
#           ground-truth cell-type labels that every method in stage 02 reads.
#
# Two kinds of script appear here:
#   .../01_preprocess/...        dataset-level preprocessing (run first)
#   .../02_methods/<M>/...       method-specific input reshaping (gene-activity
#                                matrices, region2gene links, scJoint/scBridge h5).
#                                These are preprocessing, but they only matter for
#                                the one method whose folder they sit in.
# Environment: R/Seurat/Signac for the .R steps, a general scanpy env for the .py steps.
#              See envs/ .

if want pbmc3k; then
  say "01 preprocess :: PBMC 3k  (10x multiome; 23 methods; Fig2, Fig3, Fig4, FigS1-S3)"
  note "dataset-level preprocessing"
  sub "datasets/pbmc3k/01_preprocess"
  run "datasets/pbmc3k/01_preprocess/0_pbmc3k_data_design.R"
  run "datasets/pbmc3k/01_preprocess/data_prep.R"

  note "method-specific input preparation (feeds 02_run_methods.sh)"
  sub "datasets/pbmc3k/02_methods/MaxFuse"
  run "datasets/pbmc3k/02_methods/MaxFuse/prep_maxfuse_input.R"
  sub "datasets/pbmc3k/02_methods/MinNet"
  run "datasets/pbmc3k/02_methods/MinNet/1_2_test_dataprocess.py"
  run "datasets/pbmc3k/02_methods/MinNet/1_test_dataprocess.R"
  sub "datasets/pbmc3k/02_methods/MinNet/train"
  run "datasets/pbmc3k/02_methods/MinNet/train/0_lognorm.R"
  run "datasets/pbmc3k/02_methods/MinNet/train/1_traindata_process.py"
  sub "datasets/pbmc3k/02_methods/MultiMAP"
  run "datasets/pbmc3k/02_methods/MultiMAP/0_make_atacgene.R"
  sub "datasets/pbmc3k/02_methods/scBridge"
  run "datasets/pbmc3k/02_methods/scBridge/0_scBridge_make_commengneh5.R"
  sub "datasets/pbmc3k/02_methods/scDART"
  run "datasets/pbmc3k/02_methods/scDART/make_region2gene.R"
  sub "datasets/pbmc3k/02_methods/scJoint"
  run "datasets/pbmc3k/02_methods/scJoint/0_make_scjoint_h5_control.R"
  run "datasets/pbmc3k/02_methods/scJoint/pbmc3k_scjoint_process.py"
  sub "datasets/pbmc3k/05_geneactivity_archr/01_archr"
  run "datasets/pbmc3k/05_geneactivity_archr/01_archr/archr_GeneActivity.R"
  run "datasets/pbmc3k/05_geneactivity_archr/01_archr/archr_genescores_pbmc3k.R"
fi

if want pbmc10k; then
  say "01 preprocess :: PBMC 10k (10x multiome; FigS2)"
  note "dataset-level preprocessing"
  sub "datasets/pbmc10k/01_preprocess"
  run "datasets/pbmc10k/01_preprocess/00_data_prep.R"
fi

if want pbmc_parse; then
  say "01 preprocess :: Parse PBMC (Evercode WT Mini v3; Fig6, FigS12)"
  note "dataset-level preprocessing"
  sub "datasets/pbmc_parse/01_preprocess"
  run "datasets/pbmc_parse/01_preprocess/00_data_prep.R"
  run "datasets/pbmc_parse/01_preprocess/01_azimuth_annotate_parse.R"
  run "datasets/pbmc_parse/01_preprocess/02_compare_donor_composition.R"
  run "datasets/pbmc_parse/01_preprocess/03_qc_confidence_donors.R"
  run "datasets/pbmc_parse/01_preprocess/04_build_integration_dataset.R"
fi

if want bmmc_d1; then
  say "01 preprocess :: BMMC s1d1 (public multiome; Fig6, TableS4, FigS11)"
  note "dataset-level preprocessing"
  sub "datasets/bmmc_d1/01_preprocess"
  run "datasets/bmmc_d1/01_preprocess/00_file_prep.R"

  note "barcode harmonisation; lives in 03_metrics/ but must run before the metric chain"
  sub "datasets/bmmc_d1/03_metrics"
  run "datasets/bmmc_d1/03_metrics/00_fix_bmmc_barcodes.py"
fi

if want brca; then
  say "01 preprocess :: BRCA HT243B1-S1H4 (dbGaP phs002371.v3.p1, CONTROLLED ACCESS; Fig4, FigS4-S6)"
  note "BRCA raw data is controlled access (dbGaP phs002371.v3.p1). Nothing below runs without it."
  note "dataset-level preprocessing"
  sub "datasets/brca/01_preprocess"
  run "datasets/brca/01_preprocess/00_data_prep.R"
  run "datasets/brca/01_preprocess/01_prep_HT263_train_commonpeaks.R"
  sub "datasets/brca/05_macrophage_subsample/01_preprocess"
  run "datasets/brca/05_macrophage_subsample/01_preprocess/macro_subset_rep1.R"
  run "datasets/brca/05_macrophage_subsample/01_preprocess/macro_subset_rep2.R"
  run "datasets/brca/05_macrophage_subsample/01_preprocess/macro_subset_rep3.R"
  run "datasets/brca/05_macrophage_subsample/01_preprocess/macro_subset_rep4.R"
  run "datasets/brca/05_macrophage_subsample/01_preprocess/macro_subset_rep5.R"
fi

if want rms; then
  say "01 preprocess :: RMS Mast607A (GEO GSE30978; Fig5, FigS8-S10)"
  note "dataset-level preprocessing"
  sub "datasets/rms/01_preprocess"
  lsf "datasets/rms/01_preprocess/00_cellranger_arc_sub.sh"
  run "datasets/rms/01_preprocess/01_make_common_peak_h5.R"

  note "method-specific input preparation (feeds 02_run_methods.sh)"
  sub "datasets/rms/02_methods/MaxFuse"
  run "datasets/rms/02_methods/MaxFuse/00_prep_maxfuse_input.R"
  sub "datasets/rms/02_methods/scBridge"
  run "datasets/rms/02_methods/scBridge/00_make_scbridge_h5ad.R"
  sub "datasets/rms/02_methods/scDART"
  run "datasets/rms/02_methods/scDART/00_make_region2gene.R"
  sub "datasets/rms/02_methods/scJoint"
  run "datasets/rms/02_methods/scJoint/00_make_scjoint_input.R"
fi

say "stage 1 listing complete (DRY_RUN=${DRY_RUN})"
