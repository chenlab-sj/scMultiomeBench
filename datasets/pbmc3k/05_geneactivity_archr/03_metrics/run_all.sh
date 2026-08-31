#!/bin/bash
# Score the 6 ArchR-gene-activity variants (R2.1) -> then compare vs the fig2b Signac baseline.
# PREREQ: the 6 variant latents exist (run ../run_archr_variants.sh first; compute SKIPs missing).
cd "$(dirname "$0")" || exit 1
id() { awk -F'[<>]' '/^Job </{print $2;exit}'; }
J1=$(bsub < 00_compute_metrics_sub.sh | id);                  echo "1 metrics -> $J1"
J2=$(bsub -w "done($J1)" < 02_adjust_sub.sh | id);            echo "2 adjust  -> $J2"
J3=$(bsub -w "done($J1)" < 02_peak_similarity_sub.sh | id);              echo "3 peak    -> $J3"
J4=$(bsub -w "done($J2) && done($J3)" < plot_metrics_matrix_sub.sh | id); echo "4 plot    -> $J4"
J5=$(bsub -w "done($J4)" < ../compare_sub.sh | id);        echo "5 compare -> $J5  (-> compare_geneactivity.csv + .png)"
