#!/bin/bash
# =====================================================================
# truncation.slurm -- the truncation experiment, UBC ARC (Sockeye).
# One array task = one (cell x n x rep-chunk) -> one RDS in
# simulation/output/truncation/. Checkpointed: an existing output file is
# skipped, so a re-submit resumes in seconds.
#
# WHY THIS EXISTS. The claim "truncation restores calibration but destroys
# consistency" rests on R=50, measured at POWER 0.067 (a coin flip). This
# redoes it at R=2000 (power 0.83) and additionally tests whether the result is
# TRUNCATION or KEOGH'S CAP RULE -- his quantile(w, 0.95) lands EXACTLY on the
# ~65% point mass at w==1, deleting all up-weighting while keeping all
# down-weighting. Varying `truncate` moves the cap on/off that atom.
#
#   SMOKE : sbatch --array=1-2   --export=ALL,TRUNC_R=8,TRUNC_CHUNK=8,TRUNC_B=50 truncation.slurm
#   FULL  : bash simulation/slurm/submit_truncation.sh full
#
# Array size = 3 cells x 3 n x ceil(TRUNC_R / TRUNC_CHUNK).
#   FULL: 3 x 3 x ceil(2000/32) = 9 x 63 = 567  -> --array=1-567
# =====================================================================
#SBATCH --account=YOUR_HPC_ACCOUNT
#SBATCH --job-name=steCI_trunc
#SBATCH --array=1-567
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=0:30:00
# WALLTIME: MEASURED on ARC (smoke, 32 cores, B=50): n=300 ~2 s, n=1000 ~8 s, n=5000 ~77 s per task.
# CHUNK=32 with 32 cores is ONE mclapply wave, so a 32-rep task costs about the same wall time as the
# 8-rep smoke task -- each rep gets its own core. Allowing ~3x for B=200 vs B=50, the heaviest task is
# ~2-3 min. 30 min is a ~10x margin. Do NOT inflate this: Slurm backfills SHORT jobs into gaps, so an
# over-long --time request queues far longer than it needs to (the earlier 6:00:00 was ~180x the need).
# TRUNC_CHUNK MUST stay <= --cpus-per-task. CHUNK=100 on 32 cores is exactly the bug that killed the
# main run at walltime (10 of 81 cells short, incl. three at 640 = 2x32 waves), against sim.slurm's own
# warning.
#SBATCH --output=simulation/output/logs/trunc_%A_%a.out
#SBATCH --error=simulation/output/logs/trunc_%A_%a.err
#SBATCH --mail-user=you@example.org
#SBATCH --mail-type=FAIL
# No --partition: the YOUR_HPC_ACCOUNT allocation routes to the right nodes.

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p simulation/output/truncation simulation/output/logs

module load gcc/9.4.0
module load r/4.4.0
export R_LIBS="${R_LIBS:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"

export REPO_ROOT="$SLURM_SUBMIT_DIR"
export TRUNC_OUT="$SLURM_SUBMIT_DIR/simulation/output/truncation"
export TRUNC_R="${TRUNC_R:-2000}"
export TRUNC_CHUNK="${TRUNC_CHUNK:-32}"
export TRUNC_B="${TRUNC_B:-200}"

echo "[trunc] task ${SLURM_ARRAY_TASK_ID} start $(date) on $(hostname) cores=${SLURM_CPUS_PER_TASK}"
Rscript dev/truncation_arc.R
echo "[trunc] task ${SLURM_ARRAY_TASK_ID} done  $(date)"
