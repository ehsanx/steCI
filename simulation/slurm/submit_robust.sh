#!/bin/bash
# Submit the robustness re-run array. Run from the project root ON ARC
# (never on login01/login02 -- submit only; the work happens on compute nodes).
#   bash simulation/slurm/submit_robust.sh [PROFILE]
# PROFILE in {smoke, full, misspec, tvfx, 2conf} (default: smoke).
#
# ALWAYS RUN 'smoke' FIRST. It is 9 tasks and a few minutes, and it reports the measured
# cgroup peak memory and per-replicate cost that the full profile is sized against. The
# repo's standing rule is smoke-test at R=2 before any array; this profile is that rule.
#
# NOTE ON DELIMITERS: ROBUST_STUDIES is COLON-separated. sbatch --export is itself
# comma-separated, so a comma inside a value silently truncates the list -- exactly the bug
# that made the truth run's first smoke submission a no-op for 6 of 9 tasks.
set -euo pipefail

PROFILE="${1:-smoke}"
case "$PROFILE" in
  # ---- PLUMBING ONLY. 2 replicates per cell, all three studies. Minutes.
  # Writes R2-tagged filenames, so it CANNOT overwrite production chunks.
  smoke) ARR="1-24";  TIME="0:30:00"
         EXP="ROBUST_STUDIES=2confounder:misspec:tvfx,ROBUST_REPS_OVERRIDE=2" ;;

  # ---- THE DELIVERABLE. 357 tasks: 2confounder 15, misspec 171, tvfx 171.
  full)  ARR="1-357"; TIME="0:45:00"
         EXP="ROBUST_STUDIES=2confounder:misspec:tvfx" ;;

  # ---- SINGLE STUDY. Use if you want only the cells that move visibly: misspec has the
  # largest displacement (0.79 single-draw SD at alpha0=0.40).
  misspec) ARR="1-171"; TIME="0:45:00"; EXP="ROBUST_STUDIES=misspec" ;;
  tvfx)    ARR="1-171"; TIME="0:45:00"; EXP="ROBUST_STUDIES=tvfx" ;;
  2conf)   ARR="1-15";  TIME="0:45:00"; EXP="ROBUST_STUDIES=2confounder" ;;
  *) echo "unknown profile: $PROFILE (use smoke|full|misspec|tvfx|2conf)"; exit 1 ;;
esac

mkdir -p simulation/output/logs simulation/output/robust

echo "Submitting robustness array ($PROFILE): --array=$ARR --time=$TIME --export=ALL,$EXP"
JID=$(sbatch --parsable --array="$ARR" --time="$TIME" --export=ALL,"$EXP" simulation/slurm/robust.slurm)
echo "  job id: $JID"
echo
echo "Monitor:    squeue -u \$USER"
echo "Efficiency: seff $JID"
echo "Failures:   find simulation/output/logs -name 'robust_${JID}_*.err' -size +0"
echo
echo "Checkpointed: re-running this script skips completed chunks in seconds."

case "$PROFILE" in
  smoke) cat <<'GATE'

READ THESE TWO NUMBERS BEFORE SUBMITTING 'full':

    grep -h "peak" simulation/output/logs/robust_*.out    # measured RSS -- cut --mem if small
    grep -h "wrote" simulation/output/logs/robust_*.out   # seconds per chunk -> walltime sanity

  The n=5000 chunks are 10 replicates each and are the binding constraint.
  MEASURED 2026-07-19: a 2-replicate smoke chunk at n=5000 took 129 s (64.5 s/replicate),
  so a production chunk is ~10.8 min against a 0:45:00 walltime -- 4x headroom. Peak RSS
  was 455 MB against --mem=2G.
  If a future smoke chunk at n=5000 exceeds ~4 min, the 10-replicate chunks no longer fit
  comfortably: lower chunk_reps() in dev/robust_arc.R rather than raising --time, because
  an over-long request queues far longer than it runs.

  Then DELETE the smoke output before the full run -- it is tagged R2 so it cannot
  collide, but leaving it makes the collector's completeness check noisier:

    rm -f simulation/output/robust/rob_*_R2_*.rds
GATE
  ;;
  full|misspec|tvfx|2conf) cat <<'GATE'

When it finishes, bring the per-replicate chunks home and score them THERE:

    tar czf robust_out.tgz simulation/output/robust
    # then, on the workstation, from the repo root:
    Rscript dev/robust_collect.R

  The cluster deliberately computes no truth. dev/robust_collect.R applies the exact
  truths, prints frozen-vs-exact coverage side by side, and saves the per-replicate
  matrix -- the artifact whose absence forced this whole re-run.
GATE
  ;;
esac
