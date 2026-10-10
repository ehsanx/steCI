#!/bin/bash
# Submit the E1-target (gamma_true) re-estimation array. Run from the project root
# ON ARC (never on login01/login02 -- submit only; the work happens on compute nodes).
#   bash simulation/slurm/submit_truth.sh [PROFILE]
# PROFILE in {quick, scale, smoke, full} (default: quick).
#
# MEASURED ON ARC 2026-07-19 (smoke, job 12315011/12315021, n_hd=3e4, 2 draws on 2 cores):
#   371.9 s / 551.7 s / 696.5 s at alpha.0 = 0.15 / 0.25 / 0.4  -- that is ONE wave, so it is the
#   per-draw wall cost. The same draw costs 88 / 130 / 165 s on the author's laptop, so ARC cores
#   are roughly 4x SLOWER per core. All walltimes below are sized off the ARC numbers, not the
#   laptop ones. Peak RSS was 628-985 MB per worker, i.e. ~20-32 GB at 32 workers: --mem=64G is
#   comfortable at n_hd = 3e4.
#
# NOTE ON DELIMITERS: TRUTH_NHD is COLON-separated. sbatch --export is itself comma-separated, so a
# comma inside a value silently truncates the list -- that is exactly what happened on the first
# smoke run, where the n_hd sweep never ran and 6 of 9 tasks were no-ops.
set -euo pipefail

PROFILE="${1:-quick}"
case "$PROFILE" in
  # ---- THE MAIN DELIVERABLE. Fixes the variance problem, which is what M8 is actually about.
  # 3 alpha.0 x 1 n_hd x ceil(2000/32) = 189 tasks; one mclapply wave of 32 per task.
  # WALLTIME 2:00:00, not 1:00:00. A wave of 32 was measured only at TWO workers (373/552/695 s);
  # at 32 workers these fits contend for memory bandwidth, and the same code on a laptop lost a
  # factor of ~3.4 per worker going from 1 to 10. If that carries over, a wave is ~22 min; if it is
  # worse, 1:00:00 starts killing tasks. Checkpointing makes a walltime kill recoverable but it
  # still burns the queue slot, and 189 of them is not a cheap way to learn this.
  quick) ARR="1-189"; TIME="2:00:00"
         EXP="TRUTH_K=2000,TRUTH_CHUNK=32,TRUTH_NHD=30000" ;;

  # ---- SCALING PROBE. Only needed if the n_hd sweep is going to be run: measures how aalen's cost
  # and memory actually grow with n_hd, so 'full' can be sized instead of guessed. Deliberately
  # stops at 6e4 -- at ~quadratic cost and ARC's speed, 1e5 would already be ~2 h per task.
  # 3 alpha.0 x 3 n_hd x 1 chunk = 9 tasks, 2 workers each.
  scale) ARR="1-9";   TIME="2:00:00"
         EXP="TRUTH_K=2,TRUTH_CHUNK=2,TRUTH_CORES=2,TRUTH_NHD=30000:45000:60000" ;;

  # ---- PLUMBING ONLY. 3 tasks, one n_hd, tiny. Use when changing the scripts.
  smoke) ARR="1-3";   TIME="1:00:00"
         EXP="TRUTH_K=2,TRUTH_CHUNK=2,TRUTH_CORES=2,TRUTH_NHD=30000" ;;

  # ---- THE SWEEP. Do NOT submit this until 'scale' has been read: at ARC's measured speed the
  # 3e5 arm is plausibly ~19 h per task, far beyond any sane walltime. Size it from 'scale' first,
  # and expect to reduce the top of the sweep.
  full)  ARR="1-567"; TIME="3:00:00"
         EXP="TRUTH_K=2000,TRUTH_CHUNK=32,TRUTH_NHD=30000:100000:300000" ;;
  *) echo "unknown profile: $PROFILE (use quick|scale|smoke|full)"; exit 1 ;;
esac

mkdir -p simulation/output/logs simulation/output/truth

if [ "$PROFILE" = "full" ]; then
  echo "WARNING: 'full' includes n_hd = 3e5. At ARC's measured speed that arm may need ~19 h per"
  echo "         task and far more than --mem=64G. Run 'scale' first and size it. Ctrl-C to abort."
  sleep 5
fi

echo "Submitting truth array ($PROFILE): --array=$ARR --time=$TIME --export=ALL,$EXP"
JID=$(sbatch --parsable --array="$ARR" --time="$TIME" --export=ALL,"$EXP" simulation/slurm/truth.slurm)
echo "  job id: $JID"
echo
echo "Monitor:    squeue -u \$USER"
echo "Efficiency: seff $JID"
echo "Collect:    Rscript dev/truth_collect.R     (after the array finishes)"
echo
echo "Checkpointed: re-running this script skips completed chunks in seconds."

case "$PROFILE" in
  quick) cat <<'GATE'

This is the one that matters. When it finishes:

    Rscript dev/truth_collect.R

  gives the averaged gamma_true per alpha.0 with its Monte-Carlo SE. Expect
  MCSE around 0.00014, against an E1 bias of ~0.0057 -- i.e. ~2.5%, where the
  shipped single draw sat at ~110%.
GATE
  ;;
  scale) cat <<'GATE'

When it finishes, read the cost and memory scaling:

    grep -E "MEMORY|wrote" simulation/output/logs/truth_*.out

  Fit cost vs n_hd across the three points to get the exponent, then size the
  sweep. Check peak/worker at n_hd=6e4 against --mem=64G at 32 workers.
GATE
  ;;
esac
