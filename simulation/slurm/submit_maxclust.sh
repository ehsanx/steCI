#!/bin/bash
# Submit the max.clust sensitivity array. Run from the project root ON ARC (submit only; the work
# happens on compute nodes).
#   bash simulation/slurm/submit_maxclust.sh [PROFILE]
# PROFILE in {smoke, full} (default: smoke).
#
# ALWAYS RUN 'smoke' FIRST. 4 tasks, a few minutes, and it is the only thing that MEASURES the
# seconds/fit and peak memory that FULL's walltime is sized against. Repo standing rule: smoke
# before any array.
#
# DELIMITERS: MC_CELLS is COLON-separated. sbatch --export is comma-separated, so a comma inside a
# value silently truncates the list (the bug that no-op'd the truth run's first smoke).
set -euo pipefail

PROFILE="${1:-smoke}"
case "$PROFILE" in
  # ---- PLUMBING ONLY. 2 replicates per cell at 4 cells. R2-tagged filenames cannot collide with
  # production. Worst task = 2 reps x 2 fits at n=5000 (~40 s MEASURED-ish; confirm from the log).
  smoke) ARR="1-4";   TIME="0:20:00"; MEM="4G"
         EXP="MC_CELLS=1:2:3:4,MC_REPS_OVERRIDE=2,MC_CHUNK=2" ;;

  # ---- THE DELIVERABLE. 27 cells x 100 reps, 5 reps/chunk = 540 tasks at n=5000 ONLY.
  # n<5000 is exact by construction (max.clust=1000 is not < 300/1000), so it is not run.
  # 5 reps/chunk (10 fits/task) keeps tasks SHORT: the smoke measured 29-47 s/fit on the FAST
  # evhigh cells (peak RSS ~1.3 GB); the evlow cells carry ~1.8x the stacked rows and run slower,
  # so 10 fits x ~85 s ~ 14 min sits well inside 0:45:00, and --mem=6G covers the larger evlow data.
  full)  ARR="1-540"; TIME="0:45:00"; MEM="6G"
         EXP="MC_REPS=100,MC_CHUNK=5" ;;

  *) echo "unknown profile: $PROFILE (use smoke|full)"; exit 1 ;;
esac

mkdir -p simulation/output/logs simulation/output/maxclust

echo "Submitting maxclust array ($PROFILE): --array=$ARR --time=$TIME --mem=$MEM --export=ALL,$EXP"
JID=$(sbatch --parsable --array="$ARR" --time="$TIME" --mem="$MEM" --export=ALL,"$EXP" simulation/slurm/maxclust.slurm)
echo "  job id: $JID"
echo
echo "Monitor:    squeue -u \$USER"
echo "Efficiency: seff $JID"
echo "Failures:   find simulation/output/logs -name 'mc_${JID}_*.err' -size +0"
echo
echo "Checkpointed: re-running this script skips completed chunks in seconds."

case "$PROFILE" in
  smoke) cat <<'GATE'

READ THESE TWO NUMBERS BEFORE SUBMITTING 'full':

    grep -h "peak"  simulation/output/logs/mc_*.out    # measured RSS  -- cut/raise --mem
    grep -h "wrote" simulation/output/logs/mc_*.out     # seconds/fit  -> walltime sanity

  Two aalen fits per replicate at n=5000 (the exact max.clust=NULL fit is the slower one). MEASURED
  on the evhigh smoke: 29-47 s/fit, peak RSS ~1.3 GB. The full profile uses 5 reps/chunk (10
  fits/task) at 0:45:00 / 6G because the evlow cells (~1.8x the rows) run slower and larger.

  Then DELETE the smoke output before the full run (R2-tagged, so it cannot collide, but it makes
  the collector's completeness check noisier):

    rm -f simulation/output/maxclust/mc_*_R2_*.rds
GATE
  ;;
  full) cat <<'GATE'

When it finishes, bring the per-replicate chunks home and score them THERE:

    tar czf maxclust_out.tgz simulation/output/maxclust
    # then, on the workstation, from the repo root:
    Rscript dev/maxclust_collect.R

  The cluster computes no truth. dev/maxclust_collect.R applies the exact E1 truths, prints the
  mean relative SE shift with its t-test, the per-cell shifts, and the induced coverage change,
  and writes dev/maxclust_results.rds.
GATE
  ;;
esac
