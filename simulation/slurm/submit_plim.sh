#!/bin/bash
# Submit the pooled-plim array. Run from the project root ON ARC (submit only; work happens on
# compute nodes).
#   bash simulation/slurm/submit_plim.sh [PROFILE]
# PROFILE in {smoke, full, worst} (default: smoke).
#
# ALWAYS RUN 'smoke' FIRST. 12 tasks, a few minutes, and it is the ONLY thing that measures the
# n=1e5 peak memory and seconds/fit that the full profile is sized against (the walltime and
# --mem in plim.slurm are extrapolated, not measured). Repo standing rule: smoke before any array.
#
# DELIMITERS: PLIM_CELLS / PLIM_NS are COLON-separated. sbatch --export is comma-separated, so a
# comma inside a value silently truncates the list.
set -euo pipefail

PROFILE="${1:-smoke}"
case "$PROFILE" in
  # ---- PLUMBING ONLY. 2 replicates per cell at n=1e4. 4 tasks (4 cells x 1 chunk). R2-tagged
  # filenames cannot collide with production. Worst task = 2 fits at n=1e4 (~3-4 min MEASURED).
  smoke) ARR="1-4";   TIME="0:30:00"; MEM="4G"
         EXP="PLIM_CELLS=1:2:3:4,PLIM_NS=10000,PLIM_REPS_OVERRIDE=2" ;;

  # ---- THE DELIVERABLE. 160 tasks: 4 cells x 40 chunks (5 reps/chunk) at n=1e4 ONLY.
  # n=3e4/1e5 dropped: no resolvable signal + O(n^2) cost (see dev/plim_arc.R for the
  # sizing rationale). ~22 CPU-hr total, ~10 min/task.
  full)  ARR="1-160"; TIME="0:30:00"; MEM="4G"
         EXP="PLIM_CELLS=1:2:3:4,PLIM_NS=10000" ;;

  # ---- WORST CELL ONLY (evlow_cfstrong_pvhigh), n=1e4. 40 tasks. Use to iterate cheaply.
  worst) ARR="1-40";  TIME="0:30:00"; MEM="4G"
         EXP="PLIM_CELLS=1,PLIM_NS=10000" ;;

  *) echo "unknown profile: $PROFILE (use smoke|full|worst)"; exit 1 ;;
esac

mkdir -p simulation/output/logs simulation/output/plim

echo "Submitting plim array ($PROFILE): --array=$ARR --time=$TIME --mem=$MEM --export=ALL,$EXP"
JID=$(sbatch --parsable --array="$ARR" --time="$TIME" --mem="$MEM" --export=ALL,"$EXP" simulation/slurm/plim.slurm)
echo "  job id: $JID"
echo
echo "Monitor:    squeue -u \$USER"
echo "Efficiency: seff $JID"
echo "Failures:   find simulation/output/logs -name 'plim_${JID}_*.err' -size +0"
echo
echo "Checkpointed: re-running this script skips completed chunks in seconds."

case "$PROFILE" in
  smoke) cat <<'GATE'

READ THESE TWO NUMBERS BEFORE SUBMITTING 'full':

    grep -h "peak"  simulation/output/logs/plim_*.out    # measured RSS -- cut/raise --mem
    grep -h "wrote" simulation/output/logs/plim_*.out     # seconds/fit -> walltime sanity

  n=1e4 fits are 85-122 s (MEASURED 2026-07-22), peak RSS 470-631 MB. A 5-fit task is ~10 min.

  *** n>=3e4 IS GONE ON PURPOSE. *** aalen is ~O(n^2) here (985 s/fit @3e4, ~3 h @1e5) AND the
  n>=3e4 arm carries no resolvable signal (bias moves ~0.002 vs SE ~0.008). This array is a
  finite-n=1e4 independent-code check only; the n -> infinity value is dev/plim_closed_form.R.
  Do not re-add higher n.

  Then DELETE the smoke output before the full run (it is R2-tagged so it cannot collide, but it
  makes the collector's completeness check noisier):

    rm -f simulation/output/plim/plim_*_R2_*.rds
GATE
  ;;
  full|worst) cat <<'GATE'

When it finishes, bring the per-replicate chunks home and score them THERE:

    tar czf plim_out.tgz simulation/output/plim
    # then, on the workstation, from the repo root:
    Rscript dev/plim_collect.R

  The cluster computes no truth. dev/plim_collect.R applies the exact closed-form truths, prints
  the tau=5 bias vs n (vanish-vs-plateau), and cross-checks against dev/plim_closed_form.R.
GATE
  ;;
esac
