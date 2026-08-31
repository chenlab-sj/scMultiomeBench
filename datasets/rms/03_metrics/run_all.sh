#!/bin/bash
# Full RMS Fig5a chain (dependency-chained bsub). Run on the cluster login node:
#     bash run_all.sh
#   1 compute ── 2 peak ── 3 plot
# compute scores the 14 staged rep1 latents (+ writes knn_pred_label__*.csv + adjusts accuracy); peak
# computes peakdist_adj from those predictions on the hg38 Mast607A ATAC (slow: cached bigWig export);
# plot renders the Fig5a composite table with SPLICE_PUBLISHED=0 (scores all 14 from THIS run).
cd "$(dirname "$0")" || exit 1
id(){ awk -F'[<>]' '/^Job </{print $2; exit}'; }

J1=$(bsub < 00_compute_metrics_sub.sh | id);            echo "1 compute -> $J1"
J2=$(bsub -w "done($J1)" < 02_peak_similarity_sub.sh | id);        echo "2 peak    -> $J2  (after $J1; bigWig export = slow, cached)"
J3=$(bsub -w "done($J2)" < plot_metrics_matrix_sub.sh | id);        echo "3 plot    -> $J3  (after $J2; SPLICE=0 -> Fig5a table)"

echo ""
echo "submitted. when J$J3 finishes -> metrics_matrix.csv + sum_metrics_clean.csv + metrics_table.html"
echo "  (peakdist.csv after J$J2; check figS4a_metrics.log for any SKIP/FAILED method)"
