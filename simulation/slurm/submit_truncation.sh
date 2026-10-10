#!/bin/bash
# Submit the truncation experiment array. Run from the project root ON ARC
# (never on login01/login02 -- submit only; the work happens on compute nodes).
#   bash simulation/slurm/submit_truncation.sh [PROFILE]
# PROFILE in {smoke, full} (default: smoke).
set -euo pipefail

PROFILE="${1:-smoke}"
case "$PROFILE" in
  # Array = 3 cells x 3 n x ceil(R/CHUNK). CHUNK must be <= cpus-per-task (32):
  # CHUNK=100 on 32 cores is the bug that cost the main run 10 cells at walltime.
  smoke) ARR="1-6";   EXP="TRUNC_R=8,TRUNC_CHUNK=8,TRUNC_B=50" ;;     # 3x3xceil(8/8)=9; 6 probes it
  full)  ARR="1-567"; EXP="TRUNC_R=2000,TRUNC_CHUNK=32,TRUNC_B=200" ;; # 3x3xceil(2000/32)=567
  *) echo "unknown profile: $PROFILE (use smoke|full)"; exit 1 ;;
esac

mkdir -p simulation/output/logs simulation/output/truncation

echo "Submitting truncation array ($PROFILE): --array=$ARR --export=ALL,$EXP"
JID=$(sbatch --parsable --array="$ARR" --export=ALL,"$EXP" simulation/slurm/truncation.slurm)
echo "  job id: $JID"
echo
echo "Monitor:   squeue -u \$USER"
echo "Efficiency: seff $JID"
echo "Collect:   Rscript dev/truncation_collect.R    (after the array finishes)"
echo
echo "Checkpointed: re-running this script skips completed chunks in seconds."
