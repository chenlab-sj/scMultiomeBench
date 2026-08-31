#!/bin/bash
# Fire the whole Parse cross-platform metrics pipeline with LSF job dependencies, in one shot.
# Run on the login node (it only submits jobs and returns):
#     bash run_all.sh            # full pipeline incl. the optional UMAP panels
#     bash run_all.sh --no-umap  # skip the UMAP panels
#
# Dependency shape (LSF -w "done(<id>)" = start only if the predecessor SUCCEEDED):
#     1 compute_metrics ─┬─ 2 adjust ─┐
#                        ├─ 3 peak    ┴─ 4 plot_metrics   -> the metrics figure
#                        └─ (opt) umap
# (the Parse-vs-pbmc3k rank comparison is a separate follow-up analysis: ../rank_comparison/)
cd "$(dirname "$0")" || exit 1

# pull "12345" out of: Job <12345> is submitted to queue <large_mem>.
id() { awk -F'[<>]' '/^Job </{print $2; exit}'; }

J1=$(bsub < 00_compute_metrics_sub.sh | id);                       echo "1 metrics  -> $J1"
J2=$(bsub -w "done($J1)" < 02_adjust_sub.sh | id);                 echo "2 adjust   -> $J2  (after $J1)"
J3=$(bsub -w "done($J1)" < 02_peak_similarity_sub.sh | id);                   echo "3 peak     -> $J3  (after $J1)"
J4=$(bsub -w "done($J2) && done($J3)" < plot_metrics_matrix_sub.sh | id);echo "4 plot     -> $J4  (after $J2 && $J3)"

if [ "$1" != "--no-umap" ]; then
  J5=$(bsub -w "done($J1)" < umap_sub.sh | id);                 echo "  umap     -> $J5  (after $J1)"
fi

echo ""
echo "submitted. watch:  bjobs -w    |    tail -f parse_metrics.log"
echo "final figure when J$J4 finishes: metrics_matrix.csv, sum_metrics_clean.csv, metrics_table.png"
echo "then (separate analysis): cd ../rank_comparison && bsub < rank_comparison_sub.sh"
