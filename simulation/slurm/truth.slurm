#!/bin/bash
# =====================================================================
# truth.slurm -- re-estimate the E1 target `gamma_true`, UBC ARC (Sockeye).
# One array task = one (alpha.0 x n_hd x seed-chunk) -> one RDS in
# simulation/output/truth/. Checkpointed: an existing output file is skipped,
# so a re-submit resumes in seconds.
#
# WHY THIS EXISTS. The shipped E1 target is ONE frozen Monte-Carlo draw whose
# replication SD (0.00254 / 0.00382 / 0.00629 at alpha.0 = 0.15 / 0.25 / 0.4) is
# comparable to the ~0.0057 E1 bias it is used to measure, and which -- being one
# draw shared across every replicate and scenario -- is a FIXED displacement that
# does not average out.
#
# It also sweeps n_hd, which the laptop plan could not afford and which matters:
# averaging seeds converges to E[gamma_hat(n_hd)], NOT to the population gamma.
# If the aalen const() fit is biased at n_hd = 3e4, seed-averaging will never
# reveal it. The sweep is how we find out.
#
#   SMOKE : bash simulation/slurm/submit_truth.sh smoke
#   FULL  : bash simulation/slurm/submit_truth.sh full
#
# Array size = 3 alpha.0 x |TRUTH_NHD| x ceil(TRUTH_K / TRUTH_CHUNK).
#   FULL: 3 x 3 x ceil(2000/32) = 3 x 3 x 63 = 567  -> --array=1-567
# =====================================================================
#SBATCH --account=YOUR_HPC_ACCOUNT
#SBATCH --job-name=steCI_truth
#SBATCH --array=1-567
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=1:00:00
# WALLTIME is set per-profile by submit_truth.sh --time, which overrides the value below.
# MEASURED ON ARC 2026-07-19 (smoke, n_hd=3e4, 2 draws on 2 cores): 372 / 552 / 697 s at
# alpha.0 = 0.15 / 0.25 / 0.4. That is ONE mclapply wave, so it is the per-draw wall cost, and
# CHUNK=32 on 32 cores costs about the same. ARC cores are ~4x SLOWER per core than the laptop
# these were first sized on (88 / 130 / 165 s) -- the original 1:00:00 was sized off the laptop
# and was wrong for the sweep arms.
# MEASURED memory: 628-985 MB per worker at n_hd=3e4, so ~20-32 GB at 32 workers. 64G is
# comfortable there. It has NOT been measured above n_hd=3e4; run the 'scale' profile first.
# Do NOT inflate --time blindly: Slurm backfills SHORT jobs into gaps, so an over-long request
# queues far longer than it needs to. TRUTH_CHUNK MUST stay <= --cpus-per-task.
#SBATCH --output=simulation/output/logs/truth_%A_%a.out
#SBATCH --error=simulation/output/logs/truth_%A_%a.err
#SBATCH --mail-user=you@example.org
#SBATCH --mail-type=FAIL
# No --partition: the YOUR_HPC_ACCOUNT allocation routes to the right nodes.

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p simulation/output/truth simulation/output/logs

module load gcc/9.4.0
module load r/4.4.0
export R_LIBS="${R_LIBS:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"

export REPO_ROOT="$SLURM_SUBMIT_DIR"
export TRUTH_OUT="$SLURM_SUBMIT_DIR/simulation/output/truth"
export TRUTH_K="${TRUTH_K:-2000}"
export TRUTH_CHUNK="${TRUTH_CHUNK:-32}"
export TRUTH_NHD="${TRUTH_NHD:-30000}"   # COLON-separated if a list; see submit_truth.sh
# TRUTH_CORES lets a profile use fewer workers than the allocation. The smoke profile sets it to 2
# so peak RSS at n_hd=3e5 is measured without 32 concurrent 6e5-row fits discovering the limit.
export TRUTH_CORES="${TRUTH_CORES:-$SLURM_CPUS_PER_TASK}"

echo "[truth] task ${SLURM_ARRAY_TASK_ID} start $(date) on $(hostname) cores=${SLURM_CPUS_PER_TASK}"
Rscript dev/truth_arc.R
echo "[truth] task ${SLURM_ARRAY_TASK_ID} done  $(date)"
