#!/bin/bash
# Fire the BMMC cross-site / paired metrics pipeline with LSF dependencies, in one shot.
# Run from INSIDE the experiment dir on the login node (so LS_SUBCWD resolves for each job):
#     cd BMMC_d1/benchmark/crosssite && bash run_all.sh      # or .../s1d1_paired
# PREREQUISITE: the 12 method latents must exist under .../BMMC_d1/<exp>/<method>/ (they do). compute
# SKIPS any missing latent, so you can re-run as stragglers land.
#     1 compute ─┬─ 2 major ─┐
#                └─ 3 peak   ┴─ 4 plot
cd "$(dirname "$0")" || exit 1
id() { awk -F'[<>]' '/^Job </{print $2; exit}'; }

J1=$(bsub < 00_compute_metrics_sub.sh | id);                   echo "1 compute -> $J1"
J2=$(bsub -w "done($J1)" < 03_major_recompute_sub.sh | id);    echo "2 major   -> $J2  (after $J1)"
J3=$(bsub -w "done($J1)" < 02_peak_similarity_sub.sh | id);               echo "3 peak    -> $J3  (after $J1)"
J4=$(bsub -w "done($J2) && done($J3)" < plot_metrics_matrix_sub.sh | id);  echo "4 plot    -> $J4  (after $J2 && $J3)"
echo ""
echo "submitted. watch:  bjobs -w   |   tail -f cs_metrics.log (or sp_*)"
echo "final matrix when plot finishes: metrics_matrix.csv"
