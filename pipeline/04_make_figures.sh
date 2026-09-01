#!/usr/bin/env bash
#
# scMultiomeBench -- pipeline stage 4: figure generation
#
# WHAT THIS IS
#   A thin driver that records the execution ORDER of the real analysis scripts.
#   It contains no analysis logic of its own; every line below points at a script
#   that lives under datasets/ (or summary/) and does the actual work.
#
# DRY RUN IS THE DEFAULT.
#   With DRY_RUN=1 (the default) this script only PRINTS the commands it would run.
#   Nothing is executed, no file is written, no job is submitted.
#   To actually execute:   DRY_RUN=0 ./04_make_figures.sh [dataset]
#   Do not do that casually -- see pipeline/README.md for the runtime cost.
#
# USAGE
#   ./04_make_figures.sh                 # dry run, all datasets
#   ./04_make_figures.sh pbmc3k          # dry run, one dataset
#   DRY_RUN=0 ./04_make_figures.sh rms   # actually execute, one dataset
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

# Consumes: the metric tables from stage 03 (and, for two panels only, the
#           published baseline values in results/<dataset>/published_reference/).
# Produces: the manuscript figure files.
#
# This is the cheap stage. Given the shipped metric tables it runs in minutes on a
# laptop and does not need any method environment -- an R env with ggplot2/patchwork
# and a python env with matplotlib/seaborn/scanpy cover everything here.
#
# published_reference splicing: only two scripts actually splice published baseline
# values into a panel -- pbmc3k Fig2B and the BRCA FigS6 matrix. Other scripts contain a splice
# branch that is dead code (their drivers set SPLICE_PUBLISHED=0).
#
# NOTE on Fig S7 / 4C / 5C (formerly missing, recovered 2026-08-31): Fig S7 = the three
# strip scripts in figS7_multiome_vs_annotation/ (assembled manually); Fig 4C =
# plot_fig4c_confusion.py; Fig 5C = figS7_S9/plot_pileup_grid_2page.R with the 7-row
# method selection documented in FIGURES.md.
#
# NOT IN THE REPO: predicted_ataclabel_foxo1value.csv (88.7 MB) and
# predicted_ataclabel_myod1value.csv (18.1 MB) are gitignored for size and are needed
# by the RMS FigS8/S9 pileup panels. Available from the authors on request.

if want pbmc3k; then
  say "04 figures :: PBMC 3k  (10x multiome; 23 methods; Fig2, Fig3, Fig4, FigS1-S3)"
  note "fig2b/plot_metrics_matrix.R is one of the two scripts that splice published_reference baseline values"
  sub "datasets/pbmc3k/04_figures/fig2a"
  run "datasets/pbmc3k/04_figures/fig2a/make_umap_grid.py"
  lsf "datasets/pbmc3k/04_figures/fig2a/umap_grid_sub.sh"
  sub "datasets/pbmc3k/04_figures/fig2b"
  run "datasets/pbmc3k/04_figures/fig2b/plot_metrics_matrix.R"
  lsf "datasets/pbmc3k/04_figures/fig2b/plot_metrics_matrix_sub.sh"
  run "datasets/pbmc3k/04_figures/fig2b/plot_umap.py"
  lsf "datasets/pbmc3k/04_figures/fig2b/umap_sub.sh"
  sub "datasets/pbmc3k/04_figures/fig3a"
  run "datasets/pbmc3k/04_figures/fig3a/plot_fig3a.R"
  sub "datasets/pbmc3k/04_figures/fig3bcd"
  run "datasets/pbmc3k/04_figures/fig3bcd/make_celltype_dist.py"
  sub "datasets/pbmc3k/04_figures/figS1"
  run "datasets/pbmc3k/04_figures/figS1/make_figS1_label_validation.R"
  run "datasets/pbmc3k/04_figures/figS1/marker_plot_label.R"
  lsf "datasets/pbmc3k/04_figures/figS1/marker_sub.sh"
  sub "datasets/pbmc3k/05_geneactivity_archr/04_figures"
  run "datasets/pbmc3k/05_geneactivity_archr/04_figures/compare_geneactivity.py"
  lsf "datasets/pbmc3k/05_geneactivity_archr/04_figures/compare_sub.sh"
  lsf "datasets/pbmc3k/05_geneactivity_archr/04_figures/plot_metrics_matrix_sub.sh"
fi

if want pbmc10k; then
  say "04 figures :: PBMC 10k (10x multiome; FigS2)"
  sub "datasets/pbmc10k/04_figures"
  run "datasets/pbmc10k/04_figures/plot_metrics_matrix.R"
  lsf "datasets/pbmc10k/04_figures/plot_metrics_matrix_sub.sh"
fi

if want pbmc_parse; then
  say "04 figures :: Parse PBMC (Evercode WT Mini v3; Fig6, FigS12)"
  sub "datasets/pbmc_parse/04_figures"
  run "datasets/pbmc_parse/04_figures/plot_generalization_scatter.py"
  run "datasets/pbmc_parse/04_figures/plot_metrics_matrix.R"
  lsf "datasets/pbmc_parse/04_figures/plot_metrics_matrix_sub.sh"
  run "datasets/pbmc_parse/04_figures/plot_umap.py"
  lsf "datasets/pbmc_parse/04_figures/umap_sub.sh"
fi

if want bmmc_d1; then
  say "04 figures :: BMMC s1d1 (public multiome; Fig6, TableS4, FigS11)"
  sub "datasets/bmmc_d1/04_figures"
  run "datasets/bmmc_d1/04_figures/fig6_combined.py"
  run "datasets/bmmc_d1/04_figures/fig6_panels.py"
  run "datasets/bmmc_d1/04_figures/plot_metrics_matrix.R"
  run "datasets/bmmc_d1/04_figures/tableS4_batch_wilcoxon.R"
  sub "datasets/bmmc_d1/05_crosssite/04_figures"
  run "datasets/bmmc_d1/05_crosssite/04_figures/plot_metrics_matrix.R"
  lsf "datasets/bmmc_d1/05_crosssite/04_figures/plot_metrics_matrix_sub.sh"
  sub "datasets/bmmc_d1/06_s1d1_paired/04_figures"
  run "datasets/bmmc_d1/06_s1d1_paired/04_figures/figS11_table_9metric.R"
  run "datasets/bmmc_d1/06_s1d1_paired/04_figures/plot_metrics_matrix.R"
  lsf "datasets/bmmc_d1/06_s1d1_paired/04_figures/plot_metrics_matrix_sub.sh"
fi

if want brca; then
  say "04 figures :: BRCA HT243B1-S1H4 (dbGaP phs002371.v3.p1, CONTROLLED ACCESS; Fig4, FigS4-S6)"
  note "plot_metrics_matrix.R (FigS6) is the other splice site; its _sub.sh sets SPLICE_PUBLISHED=1"
  sub "datasets/brca/04_figures"
  run "datasets/brca/04_figures/plot_fig4_combined_brca_pbmc.py"
  run "datasets/brca/04_figures/plot_fig4c_confusion.py"
  run "datasets/brca/04_figures/plot_figS4_S5_confusion.py"
  run "datasets/brca/04_figures/plot_metrics_matrix.R"
  lsf "datasets/brca/04_figures/plot_metrics_matrix_sub.sh"
  sub "datasets/brca/05_macrophage_subsample/04_figures"
  run "datasets/brca/05_macrophage_subsample/04_figures/plot_fig4d_macrophage.py"
fi

if want rms; then
  say "04 figures :: RMS Mast607A (GEO GSE209784; Fig5, FigS7-S10)"
  note "figS7_S9/plot_pileup_grid_2page.R produces FigS8-S9 and needs the two gitignored predicted_ataclabel_*.csv files"
  note "figS7_multiome_vs_annotation/ draws the three FigS7 strips (MYOD1/FOXO1/MEOX2), assembled manually"
  sub "datasets/rms/04_figures/fig5a"
  run "datasets/rms/04_figures/fig5a/plot_metrics_matrix.R"
  lsf "datasets/rms/04_figures/fig5a/plot_metrics_matrix_sub.sh"
  sub "datasets/rms/04_figures/fig5b"
  run "datasets/rms/04_figures/fig5b/plot_fig5b.py"
  sub "datasets/rms/04_figures/figS10"
  run "datasets/rms/04_figures/figS10/plot_figS10_combined.py"
  sub "datasets/rms/04_figures/figS7_S9"
  run "datasets/rms/04_figures/figS7_S9/plot_pileup_grid_2page.R"
  sub "datasets/rms/04_figures/figS7_multiome_vs_annotation"
  run "datasets/rms/04_figures/figS7_multiome_vs_annotation/plot_predpeak_final_v2.R"
  run "datasets/rms/04_figures/figS7_multiome_vs_annotation/plot_predpeak_final_foxo1_v2.R"
  run "datasets/rms/04_figures/figS7_multiome_vs_annotation/plot_predpeak_final_meox2_v2.R"
  sub "datasets/rms/04_figures/label_validation"
  run "datasets/rms/04_figures/label_validation/plot_marker_label_validation.R"
fi

if want summary; then
  say "04 figures :: cross-dataset summary (Fig7, FigS2A, FigS13)"
  note "these read metric tables from every dataset; run last"
  sub "summary/04_figures"
  run "summary/04_figures/plot_fig7_summary.R"
  lsf "summary/04_figures/plot_fig7_summary_sub.sh"
  run "summary/04_figures/plot_figS13_scalability.R"
  run "summary/04_figures/plot_figS13_weighting_rank.R"
  run "summary/04_figures/plot_figS2a_ksens_combined.py"
fi

say "stage 4 listing complete (DRY_RUN=${DRY_RUN})"
