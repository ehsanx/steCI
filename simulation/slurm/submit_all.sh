#!/bin/bash
# Submit the simulation array then the aggregation as a dependent job.
# Usage (from the project root on ARC):
#   bash simulation/slurm/submit_all.sh [PROFILE]
# PROFILE in {smoke, pilot, full} (default: pilot).
set -euo pipefail

PROFILE="${1:-pilot}"
case "$PROFILE" in
  # Jackknife (M5) is O(n) refits/replicate -> disabled here (SIM_JACK_MAXN=0);
  # run it as a separate small-n experiment. See README / PLAN.
  # Array size = 81 scenarios x ceil(N_REPS / CHUNK). Small CHUNK keeps each
  # task short so the heavy n=5000 + M4-bootstrap cells stay under --time.
  # SIM_CHUNK MUST be <= --cpus-per-task (32). CHUNK=100 on 32 cores is the bug that cost the 2026-07
  # run 10 cells at walltime -- three of them at exactly 640 reps = 2x32 mclapply waves killed mid-third
  # -- against sim.slurm's OWN warning. CHUNK=32 = one wave per task.
  smoke) ARR="1-2";    EXP="SIM_N_REPS=4,SIM_CHUNK=4,SIM_B=50,SIM_JACK_MAXN=0" ;;
  pilot) ARR="1-324";  EXP="SIM_N_REPS=100,SIM_CHUNK=32,SIM_B=200,SIM_JACK_MAXN=0" ;;   # 81 x ceil(100/32)=4
  full)  ARR="1-2592"; EXP="SIM_N_REPS=1000,SIM_CHUNK=32,SIM_B=200,SIM_JACK_MAXN=0" ;;  # 81 x ceil(1000/32)=32
  *) echo "unknown profile: $PROFILE (use smoke|pilot|full)"; exit 1 ;;
esac

# Slurm needs the --output/--error directories to exist at submission time.
mkdir -p simulation/output/logs simulation/output/intermediate simulation/output/results

echo "Submitting sim array ($PROFILE): --array=$ARR --export=ALL,$EXP"
SIM_JID=$(sbatch --parsable --array="$ARR" --export=ALL,"$EXP" simulation/slurm/sim.slurm)
echo "  sim job id: $SIM_JID"

# afterany (not afterok): aggregate runs once the array *finishes*, even if some
# tasks failed -- aggregate.R tolerates partial/missing cells. Avoids the
# DependencyNeverSatisfied trap when a single task times out.
AGG_JID=$(sbatch --parsable --dependency=afterany:"$SIM_JID" simulation/slurm/aggregate.slurm)
echo "  aggregate job id: $AGG_JID (runs after $SIM_JID finishes, pass or fail)"
echo "Monitor with:  squeue -u \$USER   |   seff $SIM_JID"
