#!/usr/bin/env bash
#
# multiomeBench -- pipeline stage 3: metric computation
#
# WHAT THIS IS
#   A thin driver that records the execution ORDER of the real analysis scripts.
#   It contains no analysis logic of its own; every line below points at a script
#   that lives under datasets/ (or summary/) and does the actual work.
#
# DRY RUN IS THE DEFAULT.
#   With DRY_RUN=1 (the default) this script only PRINTS the commands it would run.
#   Nothing is executed, no file is written, no job is submitted.
#   To actually execute:   DRY_RUN=0 ./03_compute_metrics.sh [dataset]
#   Do not do that casually -- see pipeline/README.md for the runtime cost.
#
# USAGE
#   ./03_compute_metrics.sh                 # dry run, all datasets
#   ./03_compute_metrics.sh pbmc3k          # dry run, one dataset
#   DRY_RUN=0 ./03_compute_metrics.sh rms   # actually execute, one dataset
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

# Consumes: the latent embeddings written by stage 02, plus the ground-truth
#           labels from stage 01.
# Produces: the per-dataset metric tables that stage 04 turns into figures.
#
# Order matters and is encoded in the 00_/01_/02_... filename prefixes; the
# listing below is sorted by path so those prefixes read in order.
#
# Composite score: 6 metrics on unpaired data
#   ks_celltype_mean, asw, omics_asw, ks_inter_celltype_mean, average_accu, peakdist_adj
# Paired data adds ks.statistic, ari, ami -> 9 metrics. ARI/AMI and the omics KS
# are same-cell metrics, so they are undefined for independently sequenced datasets.
#
# KNN label transfer: k = 10, cosine distance, distance-weighted.
#
# Environments: the .py metric steps run in a general scanpy/scikit-learn env;
# the *02_peak_similarity.R steps need R with the Signac/GenomicRanges stack and read
# common/plot_regionpeak_fun.R via $BENCHMARK_FUN_DIR. NOTE: common/ is empty in the
# current tree -- the three shared helpers (benchmark_fun.py, export_groupbwg.R,
# plot_regionpeak_fun.R) are added by a later stage of the repository build.

if want pbmc3k; then
  say "03 metrics :: PBMC 3k  (10x multiome; 23 methods; Fig2, Fig3, Fig4, FigS1-S3)"
  sub "datasets/pbmc3k/03_metrics"
  run "datasets/pbmc3k/03_metrics/00_compute_metrics.py"
  lsf "datasets/pbmc3k/03_metrics/00_compute_metrics_sub.sh"
  run "datasets/pbmc3k/03_metrics/01_adjust_accuracy.py"
  run "datasets/pbmc3k/03_metrics/02_peak_similarity.R"
  lsf "datasets/pbmc3k/03_metrics/02_peak_similarity_sub.sh"
  sub "datasets/pbmc3k/03_metrics/ksensitivity"
  run "datasets/pbmc3k/03_metrics/ksensitivity/00_compute_ktest_legacy18.py"
  lsf "datasets/pbmc3k/03_metrics/ksensitivity/00_compute_ktest_legacy18_sub.sh"
  run "datasets/pbmc3k/03_metrics/ksensitivity/01_compute_ktest_new_methods.py"
  lsf "datasets/pbmc3k/03_metrics/ksensitivity/01_compute_ktest_new_methods_sub.sh"
  sub "datasets/pbmc3k/03_metrics/label_sensitivity"
  run "datasets/pbmc3k/03_metrics/label_sensitivity/00_singleR_annotate.R"
  lsf "datasets/pbmc3k/03_metrics/label_sensitivity/00_singleR_annotate_sub.sh"
  run "datasets/pbmc3k/03_metrics/label_sensitivity/01_singleR_confidence.R"
  run "datasets/pbmc3k/03_metrics/label_sensitivity/02_remap_labels.R"
  lib "datasets/pbmc3k/03_metrics/label_sensitivity/fine_label_map.R"
  sub "datasets/pbmc3k/03_metrics/label_sensitivity/rank_stability"
  run "datasets/pbmc3k/03_metrics/label_sensitivity/rank_stability/rank_stability.py"
  lsf "datasets/pbmc3k/03_metrics/label_sensitivity/rank_stability/recompute_sub.sh"
  sub "datasets/pbmc3k/03_metrics/reproducibility"
  run "datasets/pbmc3k/03_metrics/reproducibility/00_reproduce_curated.py"
  lsf "datasets/pbmc3k/03_metrics/reproducibility/00_reproduce_curated_sub.sh"
  sub "datasets/pbmc3k/05_geneactivity_archr/03_metrics"
  run "datasets/pbmc3k/05_geneactivity_archr/03_metrics/00_fix_latent_format.py"
  lsf "datasets/pbmc3k/05_geneactivity_archr/03_metrics/01_compute_metrics_sub.sh"
  lsf "datasets/pbmc3k/05_geneactivity_archr/03_metrics/02_adjust_sub.sh"
  lsf "datasets/pbmc3k/05_geneactivity_archr/03_metrics/03_peak_sub.sh"
  lsf "datasets/pbmc3k/05_geneactivity_archr/03_metrics/run_all.sh"
fi

if want pbmc10k; then
  say "03 metrics :: PBMC 10k (10x multiome; FigS2)"
  sub "datasets/pbmc10k/03_metrics"
  run "datasets/pbmc10k/03_metrics/00_compute_metrics.py"
  lsf "datasets/pbmc10k/03_metrics/00_compute_metrics_sub.sh"
  run "datasets/pbmc10k/03_metrics/01_adjust_accuracy.py"
  run "datasets/pbmc10k/03_metrics/02_peak_similarity.R"
  lsf "datasets/pbmc10k/03_metrics/02_peak_similarity_sub.sh"
  run "datasets/pbmc10k/03_metrics/03_compute_ktest_new_methods.py"
  lsf "datasets/pbmc10k/03_metrics/03_compute_ktest_new_methods_sub.sh"
  run "datasets/pbmc10k/03_metrics/04_compute_ktest_reverse.py"
  lsf "datasets/pbmc10k/03_metrics/04_compute_ktest_reverse_sub.sh"
fi

if want pbmc_parse; then
  say "03 metrics :: Parse PBMC (Evercode WT Mini v3; Fig6, FigS12)"
  sub "datasets/pbmc_parse/03_metrics"
  run "datasets/pbmc_parse/03_metrics/00_build_label.py"
  run "datasets/pbmc_parse/03_metrics/01_prep_latents.py"
  run "datasets/pbmc_parse/03_metrics/02_compute_metrics.py"
  lsf "datasets/pbmc_parse/03_metrics/02_compute_metrics_sub.sh"
  run "datasets/pbmc_parse/03_metrics/03_adjust_accuracy.py"
  lsf "datasets/pbmc_parse/03_metrics/03_adjust_sub.sh"
  run "datasets/pbmc_parse/03_metrics/04_peak_similarity.R"
  lsf "datasets/pbmc_parse/03_metrics/04_peak_sub.sh"
  lsf "datasets/pbmc_parse/03_metrics/run_all.sh"
fi

if want bmmc_d1; then
  say "03 metrics :: BMMC s1d1 (public multiome; Fig6, TableS4, FigS11)"
  sub "datasets/bmmc_d1/03_metrics"
  run "datasets/bmmc_d1/03_metrics/01_prep_latents.py"
  lsf "datasets/bmmc_d1/03_metrics/02_compute_sub.sh"
  run "datasets/bmmc_d1/03_metrics/02_run_bmmc_metrics.py"
  run "datasets/bmmc_d1/03_metrics/03_knn_bmmc.py"
  run "datasets/bmmc_d1/03_metrics/04_adjust_accuracy.py"
  run "datasets/bmmc_d1/03_metrics/05_peak_similarity.R"
  lsf "datasets/bmmc_d1/03_metrics/05_peak_sub.sh"
  lsf "datasets/bmmc_d1/03_metrics/06_portal_append_sub.sh"
  lib "datasets/bmmc_d1/03_metrics/benchmark_metrics_lib.py"
  sub "datasets/bmmc_d1/05_crosssite/03_metrics"
  run "datasets/bmmc_d1/05_crosssite/03_metrics/00_build_label.py"
  run "datasets/bmmc_d1/05_crosssite/03_metrics/01_prep_latents.py"
  run "datasets/bmmc_d1/05_crosssite/03_metrics/02_compute_metrics.py"
  lsf "datasets/bmmc_d1/05_crosssite/03_metrics/02_compute_metrics_sub.sh"
  lsf "datasets/bmmc_d1/05_crosssite/03_metrics/03_major_recompute_sub.sh"
  run "datasets/bmmc_d1/05_crosssite/03_metrics/04_adjust_accuracy.py"
  run "datasets/bmmc_d1/05_crosssite/03_metrics/05_peak_similarity.R"
  lsf "datasets/bmmc_d1/05_crosssite/03_metrics/05_peak_sub.sh"
  lsf "datasets/bmmc_d1/05_crosssite/03_metrics/run_all.sh"
  sub "datasets/bmmc_d1/06_s1d1_paired/03_metrics"
  run "datasets/bmmc_d1/06_s1d1_paired/03_metrics/00_build_label.py"
  run "datasets/bmmc_d1/06_s1d1_paired/03_metrics/01_prep_latents.py"
  run "datasets/bmmc_d1/06_s1d1_paired/03_metrics/02_compute_metrics.py"
  lsf "datasets/bmmc_d1/06_s1d1_paired/03_metrics/02_compute_metrics_sub.sh"
  lsf "datasets/bmmc_d1/06_s1d1_paired/03_metrics/03_major_recompute_sub.sh"
  run "datasets/bmmc_d1/06_s1d1_paired/03_metrics/04_adjust_accuracy.py"
  run "datasets/bmmc_d1/06_s1d1_paired/03_metrics/05_peak_similarity.R"
  lsf "datasets/bmmc_d1/06_s1d1_paired/03_metrics/05_peak_sub.sh"
  lsf "datasets/bmmc_d1/06_s1d1_paired/03_metrics/run_all.sh"
fi

if want brca; then
  say "03 metrics :: BRCA HT243B1-S1H4 (dbGaP phs002371.v3.p1, CONTROLLED ACCESS; Fig4, FigS4-S6)"
  sub "datasets/brca/03_metrics"
  run "datasets/brca/03_metrics/00_compute_metrics.py"
  lsf "datasets/brca/03_metrics/00_compute_metrics_sub.sh"
  lsf "datasets/brca/03_metrics/01_major_recompute_sub.sh"
  run "datasets/brca/03_metrics/02_adjust_accuracy.py"
  lsf "datasets/brca/03_metrics/02_adjust_accuracy_sub.sh"
  run "datasets/brca/03_metrics/03_peak_similarity.R"
  lsf "datasets/brca/03_metrics/03_peak_similarity_sub.sh"
  sub "datasets/brca/03_metrics/reproducibility"
  run "datasets/brca/03_metrics/reproducibility/00_cobolt_nmi.py"
  lsf "datasets/brca/03_metrics/reproducibility/00_cobolt_nmi_sub.sh"
  run "datasets/brca/03_metrics/reproducibility/01_midas_pick.py"
  lsf "datasets/brca/03_metrics/reproducibility/01_midas_pick_sub.sh"
  run "datasets/brca/03_metrics/reproducibility/02_reproduce_metrics.py"
  lsf "datasets/brca/03_metrics/reproducibility/02_reproduce_metrics_sub.sh"
  run "datasets/brca/03_metrics/reproducibility/03_fig4_reproduce_7type.ipynb"
  lsf "datasets/brca/03_metrics/reproducibility/03_fig4_reproduce_nb_sub.sh"
  sub "datasets/brca/05_macrophage_subsample/03_metrics"
  run "datasets/brca/05_macrophage_subsample/03_metrics/compute_new_methods_accu.py"
fi

if want rms; then
  say "03 metrics :: RMS Mast607A (GEO GSE30978; Fig5, FigS8-S10)"
  sub "datasets/rms/03_metrics"
  run "datasets/rms/03_metrics/00_prep_latents.py"
  run "datasets/rms/03_metrics/01_compute_metrics.py"
  lsf "datasets/rms/03_metrics/01_compute_metrics_sub.sh"
  run "datasets/rms/03_metrics/02_adjust_accuracy.py"
  run "datasets/rms/03_metrics/03_peak_similarity.R"
  lsf "datasets/rms/03_metrics/03_peak_similarity_sub.sh"
  run "datasets/rms/03_metrics/04_compute_new_pileups.R"
  lsf "datasets/rms/03_metrics/04_compute_new_pileups_sub.sh"
  sub "datasets/rms/03_metrics/reproducibility"
  run "datasets/rms/03_metrics/reproducibility/00_pick3_reps.py"
  lsf "datasets/rms/03_metrics/reproducibility/00_pick3_reps_scjoint7_sub.sh"
  lsf "datasets/rms/03_metrics/reproducibility/00_pick3_reps_sub.sh"
  run "datasets/rms/03_metrics/reproducibility/01_select_pick3_triplets.py"
  run "datasets/rms/03_metrics/reproducibility/02_reproduce_metrics.py"
  lsf "datasets/rms/03_metrics/reproducibility/02_reproduce_metrics_sub.sh"
  sub "datasets/rms/03_metrics"
  lsf "datasets/rms/03_metrics/run_all.sh"
  sub "datasets/rms/05_atac_precede_rna/03_metrics"
  run "datasets/rms/05_atac_precede_rna/03_metrics/00_dpt_pseudotime.py"
  run "datasets/rms/05_atac_precede_rna/03_metrics/01_region_atac_by_bin.R"
  lsf "datasets/rms/05_atac_precede_rna/03_metrics/02_run_velocyto_sub.sh"
fi

if want summary; then
  say "03 metrics :: cross-dataset summary"
  note "run only after every per-dataset metric table above exists"
  sub "summary/03_metrics"
  run "summary/03_metrics/00_weighting_sensitivity.R"
fi

say "stage 3 listing complete (DRY_RUN=${DRY_RUN})"
