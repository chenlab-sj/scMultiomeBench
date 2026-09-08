#!/bin/bash
# Fire the whole HT137 (Fig S4B) UNPAIRED metrics pipeline with LSF job dependencies, in one shot.
# Run on the cluster login node (it only submits jobs and returns):
#     bash run_all.sh            # full pipeline incl. the optional UMAP panels
#     bash run_all.sh --no-umap  # skip the UMAP panels
#
# PREREQUISITE: the method latents must exist. The 11 carried-over methods are already present; run
# the 3 new methods + scVI reps first via ../../HT137B1-S1H7/run_methods.sh (compute SKIPS any latent
# not yet there, so you can also run this now for the 11 and re-run once the new latents land).
#
# Dependency shape (LSF -w "done(<id>)" = start only if the predecessor SUCCEEDED):
#     1 compute (build_label + prep_latents + score) ─┬─ 2 adjust ─┐
#                                                      ├─ 3 peak    ┴─ 4 plot  -> the FigS4B matrix
#                                                      └─ (opt) umap
# (the HT137-vs-same-platform rank comparison is a separate follow-up: ./rank_comparison/)
cd "$(dirname "$0")" || exit 1
id() { awk -F'[<>]' '/^Job </{print $2; exit}'; }

J1=$(bsub < compute_metrics_sub.sh | id);                          echo "1 metrics -> $J1"
J2=$(bsub -w "done($J1)" < adjust_sub.sh | id);                    echo "2 adjust  -> $J2  (after $J1)"
J3=$(bsub -w "done($J1)" < peak_sub.sh | id);                      echo "3 peak    -> $J3  (after $J1)"
J4=$(bsub -w "done($J2) && done($J3)" < plot_sub.sh | id);         echo "4 plot    -> $J4  (after $J2 && $J3)"

if [ "$1" != "--no-umap" ]; then
  J5=$(bsub -w "done($J1)" < umap_sub.sh | id);                    echo "  umap    -> $J5  (after $J1)"
fi

echo ""
echo "submitted. watch:  bjobs -w   |   tail -f figS4b_metrics.log"
echo "final figure when J$J4 finishes: metrics_matrix.csv, sum_metrics_clean.csv, metrics_table.png"
