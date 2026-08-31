#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

# Submit all 25 MaxFuse macrophage-subsample jobs (5 macrosub x 5 sub) as SEPARATE large_mem CPU jobs ->
# they run in PARALLEL (and don't compete with the dgx scButterfly jobs). Each job does the 2-stage
# pipeline (Signac prep -> Fusor) writing only to {MACROSUB}/maxfuse/sub{SUB}/. Usage: bash 02_submit_maxfuse_subset.sh
SUBS=${PROJECT_ROOT}/BRCA/HT243-S1H4_subsample
n=0
for M in HT243_S1H4_macrosub HT243_S1H4_macrosub2 HT243_S1H4_macrosub3 HT243_S1H4_macrosub4 HT243_S1H4_macrosub5; do
  for N in 1 2 3 4 5; do
    O="$SUBS/$M/maxfuse/sub$N"; mkdir -p "$O"
    bsub -P benchmark -J "mf_${M#HT243_S1H4_}_s${N}" -q large_mem -R "rusage[mem=200000]" -n 4 \
         -o "$O/run.log" -e "$O/run.err" \
         "export MACROSUB=$M SUB=$N; bash $SUBS/run_maxfuse_subset_job.sh"
    n=$((n + 1))
  done
done
echo "submitted $n MaxFuse subsample jobs to large_mem (parallel; outputs in <macrosub>/maxfuse/sub<N>/)"
