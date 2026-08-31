#!/bin/bash
# Load path configuration (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT, ...)
[ -f "${CONFIG_SH:-$(git rev-parse --show-toplevel 2>/dev/null)/config/config.sh}" ] && . "${CONFIG_SH:-$(git rev-parse --show-toplevel)/config/config.sh}"

# ============================================================================================
# Submit ALL additional reproducibility reps (rep4 + rep5) for the 5 methods, with LSF job
# dependencies enforcing the required ordering, then a final 5-rep recompute gated on all of them.
#
#   Normal run (submit reps + gated recompute):   bash submit_all_reps45.sh
#   Recompute ONLY (reps already submitted):      RECOMPUTE_ONLY=1 bash submit_all_reps45.sh
#
# The recompute is gated:
#   * normal mode  -> on the JOB IDs captured here (only the ones that actually submitted);
#   * RECOMPUTE_ONLY -> on the rep JOB NAMES from your earlier run (so it waits for jobs already
#     queued). Use this after a partial run instead of re-submitting every rep.
#
# Ordering (why deps exist): within a method the reps share fixed LSA working files, so rep5 waits
# for rep4 (scDART Examples staging / scBridge latent.csv+pkl / scJoint 02_config.py+output / Portal
# staging); MIDAS runs both reps serially in one job. Different METHODS run in parallel.
# ============================================================================================

SCRIPTS=${PROJECT_ROOT}/pbmc/pbmc3k/scripts
PORTAL=${TOOLS_ROOT}/portal/Portal/pbmc3k
FIG4=${PROJECT_ROOT}/pbmc/pbmc3k/benchmark/fig4

# submit <full_script_path> [dep_jobid] -> prints JOBID on stdout ("" on failure); log to stderr.
# If dep is given but EMPTY (upstream submit failed), submits WITHOUT -w and warns, so the chain
# still runs rather than cascading into more "No matching job found" errors.
submit () {
    local script="$1" dep="$2" dir base out jid
    dir=$(dirname "$script"); base=$(basename "$script")
    if ! cd "$dir" 2>/dev/null; then echo "  !! cd failed: $dir" >&2; echo ""; return; fi
    if [ -n "$dep" ]; then
        out=$(bsub -w "done($dep)" < "$base" 2>&1)
    else
        out=$(bsub < "$base" 2>&1)
    fi
    jid=$(printf '%s\n' "$out" | grep -oE 'Job <[0-9]+>' | grep -oE '[0-9]+' | head -1)
    if [ -z "$jid" ]; then
        echo "  !! SUBMIT FAILED: ${dir##*/pbmc3k/}/$base ${dep:+(dep $dep)} :: $(printf '%s' "$out" | tr '\n' ' ')" >&2
    else
        echo "  ${dir##*/pbmc3k/}/$base  ->  job $jid${dep:+   (after $dep)}" >&2
    fi
    echo "$jid"
}

# gated_recompute <dep_expr>  (dep_expr may be empty -> submit ungated)
gated_recompute () {
    local depexpr="$1" out jid
    cd "$FIG4" || { echo "  !! cd failed: $FIG4" >&2; return; }
    if [ -n "$depexpr" ]; then out=$(bsub -w "$depexpr" < extrareps_metrics_sub.sh 2>&1)
    else                        out=$(bsub < extrareps_metrics_sub.sh 2>&1); fi
    jid=$(printf '%s\n' "$out" | grep -oE 'Job <[0-9]+>' | grep -oE '[0-9]+' | head -1)
    [ -n "$jid" ] && echo "  recompute -> job $jid   (after: ${depexpr:-<none>})" >&2 \
                  || echo "  !! recompute submit FAILED :: $(printf '%s' "$out" | tr '\n' ' ')" >&2
}

# ---- RECOMPUTE-ONLY: reps already submitted; gate on their JOB NAMES ----------------------------
if [ -n "$RECOMPUTE_ONLY" ]; then
    echo ">>> RECOMPUTE ONLY: gating on rep job NAMES from the previous run"
    NAMEDEP='done(pbmc3k_scdart_rep5) && done(pbmc3k_scb_s2_rep5) && done(pbmc3k_scjoint_s3_rep5) && done(pbmc3k_portal_rep5) && done(pbmc3k_midas_reps45) && done(scBridge_base_lat)'
    gated_recompute "$NAMEDEP"
    echo "(if any of those jobs already finished/were renamed, drop it from NAMEDEP and re-run)"
    exit 0
fi

# ---- NORMAL: submit every rep, then recompute gated on whatever actually submitted -------------
echo ">>> scDART (rep4 -> rep5)"
sd4=$(submit "$SCRIPTS/scDART/rep4/scdart_sub.sh")
sd5=$(submit "$SCRIPTS/scDART/rep5/scdart_sub.sh" "$sd4")

echo ">>> scBridge (rep4 s1->s2 -> rep5 s1->s2)"
sb4a=$(submit "$SCRIPTS/scBridge/rep4/scBridge_s1_sub.sh")
sb4b=$(submit "$SCRIPTS/scBridge/rep4/scBridge_s2_sub.sh" "$sb4a")
sb5a=$(submit "$SCRIPTS/scBridge/rep5/scBridge_s1_sub.sh" "$sb4b")
sb5b=$(submit "$SCRIPTS/scBridge/rep5/scBridge_s2_sub.sh" "$sb5a")

echo ">>> scJoint (rep4 -> rep5)"
sj4=$(submit "$SCRIPTS/scJoint/rep4/scjoint_s3_sub.sh")
sj5=$(submit "$SCRIPTS/scJoint/rep5/scjoint_s3_sub.sh" "$sj4")

echo ">>> Portal (rep4 -> rep5)"
pt4=$(submit "$PORTAL/rep4/portal_sub.sh")
pt5=$(submit "$PORTAL/rep5/portal_sub.sh" "$pt4")

echo ">>> MIDAS (rep4+rep5 serial in one job)"
md=$(submit "$SCRIPTS/midas/reps45_sub.sh")

echo ">>> scBridge rep1 recovery (extract_base)"
eb=$(submit "$SCRIPTS/scBridge/extract_base_sub.sh")

# Build the recompute dependency from ONLY the job IDs that actually came back (skip empties).
echo ">>> gated 5-rep recompute"
depexpr=""
for pair in "scDART:$sd5" "scBridge:$sb5b" "scJoint:$sj5" "Portal:$pt5" "MIDAS:$md" "extract_base:$eb"; do
    id=${pair##*:}; name=${pair%%:*}
    if [ -n "$id" ]; then depexpr="${depexpr:+$depexpr && }done($id)"
    else echo "  !! ${name}: no job id -> recompute will NOT wait on it (submit failed above)" >&2; fi
done
gated_recompute "$depexpr"

echo ""
echo "Submitted. Watch: bjobs -w   |   results -> fig4/extrareps/reproducibility.csv + compare_3v5.csv"
