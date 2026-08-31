#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

# Submit all 25 scButterfly macrophage-subsample jobs (5 macrosub reps x 5 sub-levels) as SEPARATE dgx
# GPU jobs -> they run in PARALLEL. Each writes only to {MACROSUB}/scbutterfly/sub{SUB}/ (model, latent,
# runtime all isolated), so simultaneous runs never collide. Re-runnable: a job that already trained
# reloads its saved model and just re-embeds. Usage:  bash submit_scbutterfly_subset.sh
SUBS=${PROJECT_ROOT}/BRCA/HT243-S1H4_subsample
SCRIPT="$SUBS/run_scbutterfly_subset.py"
n=0
for M in HT243_S1H4_macrosub HT243_S1H4_macrosub2 HT243_S1H4_macrosub3 HT243_S1H4_macrosub4 HT243_S1H4_macrosub5; do
  for N in 1 2 3 4 5; do
    O="$SUBS/$M/scbutterfly/sub$N"; mkdir -p "$O"
    bsub -P benchmark -J "scb_${M#HT243_S1H4_}_s${N}" -q dgx -R "rusage[mem=200000]" -gpu "num=1" \
         -o "$O/run.log" -e "$O/run.err" \
         "module load conda3/202311; source activate scbutterfly-env; export PYTHONNOUSERSITE=1 MACROSUB=$M SUB=$N; python $SCRIPT"
    n=$((n + 1))
  done
done
echo "submitted $n scButterfly subsample jobs to dgx (parallel; outputs in <macrosub>/scbutterfly/sub<N>/)"
