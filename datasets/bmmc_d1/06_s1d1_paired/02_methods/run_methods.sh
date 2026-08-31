#!/bin/bash
# Submit all 12 BMMC integration method runs for THIS experiment (the folder this script lives in).
# Run on the cluster login node:   bash run_methods.sh
# Each job is submitted FROM its own method dir (LSF runs the job in the submission cwd, so the relative
# `python run_*.py` / `R CMD BATCH *.R` resolve). Multi-stage methods are chained with LSF done() deps.
# Outputs land in <repo>/BMMC_d1/<this-experiment>/<method>/.  Then run the metrics pipeline.
cd "$(dirname "$0")" || exit 1
id()     { awk -F'[<>]' '/^Job </{print $2; exit}'; }
sub()    { ( cd "$1" && bsub < "$2" ) | id; }                 # $1=method dir  $2=sub.sh
subdep() { ( cd "$1" && bsub -w "done($3)" < "$2" ) | id; }   # $3=predecessor job id

echo "== single-stage methods =="
for spec in \
  "scglue|scglue_sub.sh" "scglue_paired|scglue_sub.sh" "scVI|01_run_scvi_sub.sh" "midas|run_midas_sub.sh" \
  "simba|run_simba_sub.sh" "seurat3|seurat3_sub.sh" "scbutterfly|run_scbutterfly_sub.sh" \
  "portal|run_portal_sub.sh" "bindsc|bindsc_sub.sh" "maxfuse|01_run_maxfuse_sub.sh"; do
  d="${spec%%|*}"; s="${spec#*|}"
  printf "  %-14s -> %s\n" "$d" "$(sub "$d" "$s")"
done

echo "== scBridge: s0 -> s1 -> s2 =="
b0=$(sub scBridge 00_scBridge_s0_sub.sh);           echo "  s0 -> $b0"
b1=$(subdep scBridge 01_scBridge_s1_sub.sh "$b0");  echo "  s1 -> $b1 (after $b0)"
b2=$(subdep scBridge 02_scBridge_s2_sub.sh "$b1");  echo "  s2 -> $b2 (after $b1)"

echo "== scjoint: s1 -> s2 -> patch(config) -> s3(+s4 assemble) =="
j1=$(sub scjoint 00_scjoint_s1_sub.sh);              echo "  s1    -> $j1"
j2=$(subdep scjoint 01_scjoint_s2_sub.sh "$j1");     echo "  s2    -> $j2 (after $j1)"
jp=$(subdep scjoint 02_scjoint_patch_sub.sh "$j2");  echo "  patch -> $jp (after $j2)"
j3=$(subdep scjoint 02_scjoint_s3_sub.sh "$jp");     echo "  s3+s4 -> $j3 (after $jp)"

echo
echo "submitted. watch:  bjobs -w"
echo "latents land in <repo>/BMMC_d1/<experiment>/<method>/ ; then run the metrics pipeline."
I need to help me polish Table S1 , also do you think I can use Table S2 to address this comments, do you agree with my response  and do you think this is enough