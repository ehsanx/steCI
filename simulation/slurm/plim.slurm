#!/bin/bash
# =====================================================================
# plim.slurm -- empirical population limit of the POOLED STE estimator, UBC ARC (Sockeye).
# One array task = one (cell x n x replicate-chunk) -> one RDS in simulation/output/plim/.
# Checkpointed: an existing output file is skipped, so a re-submit resumes in seconds.
#
# WHY THIS EXISTS. This checks whether the residual long-horizon E2 bias is a
# NON-VANISHING pooled-origin offset or a finite-sample effect. Three sample sizes cannot tell
# a vanishing bias from a plateau, so this extends the trend to large n (up to 1e5) at high
# precision. The empirical plim is then cross-checked against the closed-form pooled plim
# (dev/plim_closed_form.R); agreement pins the pooled offset exactly.
#
# THIS JOB COMPUTES NO TRUTH. It emits per-replicate (gamma_hat, mrd_hat[1..5]) via the SHIPPED
# ste_estimate() -- the estimator is not reimplemented -- and stops. Scoring against the exact
# truths, the vanish-vs-plateau read, and the closed-form cross-check happen on the workstation
# in dev/plim_collect.R. Any future truth revision is then a seconds-long re-score.
#
#   SMOKE : bash simulation/slurm/submit_plim.sh smoke     # 12 tasks (2 reps/cell-n), minutes
#   FULL  : bash simulation/slurm/submit_plim.sh full      # 152 tasks
# =====================================================================
#SBATCH --account=YOUR_HPC_ACCOUNT
#SBATCH --job-name=steCI_plim
#SBATCH --array=1-160
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=0:45:00
# WALLTIME/ARRAY/MEM are set per-profile by submit_plim.sh, which overrides the values above.
#
# SIZING -- MEASURED ON ARC (smoke, job 12342697, 2026-07-22), not extrapolated.
#   Per-fit wall cost: 85-122 s at n=1e4. timereg::aalen is empirically ~O(n^2) here (985 s/fit at
#   n=3e4, ~3 h at n=1e5), so n IS CAPPED AT 1e4. n=3e4 and n=1e5 were dropped for TWO reasons:
#   the O(n^2) cost (all four n=1e5 smoke tasks were cancelled at the walltime), AND -- decisively --
#   the n>=3e4 arm carries no resolvable signal (bias moves ~0.002 vs an SE ~0.008; see dev/plim_arc.R).
#   This array now does one thing: an independent-code (real aalen) check
#   that the finite-n bias at n=1e4 is clearly positive. The n -> infinity value is the closed form
#   (dev/plim_closed_form.R), not this array.
#   160 tasks x 5 fits @ n=1e4 = ~10 min/task; 0:30:00 is ~3x headroom.
#
# MEMORY -- MEASURED: peak RSS 470-631 MB at n=1e4. --mem=4G is ~6x headroom.
#SBATCH --output=simulation/output/logs/plim_%A_%a.out
#SBATCH --error=simulation/output/logs/plim_%A_%a.err
#SBATCH --mail-user=you@example.org
#SBATCH --mail-type=FAIL
# No --partition: the YOUR_HPC_ACCOUNT allocation routes to the right nodes.

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p simulation/output/plim simulation/output/logs

module load gcc/9.4.0
module load r/4.4.0
export R_LIBS="${R_LIBS:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"

export REPO_ROOT="$SLURM_SUBMIT_DIR"
export PLIM_OUT="$SLURM_SUBMIT_DIR/simulation/output/plim"
# PLIM_CELLS and PLIM_NS are COLON-separated. sbatch --export is comma-separated, so a comma
# inside a value silently truncates the list -- the bug that no-op'd the truth run's first smoke.
export PLIM_CELLS="${PLIM_CELLS:-1:2:3:4}"
export PLIM_NS="${PLIM_NS:-10000:30000}"

echo "[plim] task ${SLURM_ARRAY_TASK_ID} start $(date) on $(hostname)"
Rscript dev/plim_arc.R
echo "[plim] task ${SLURM_ARRAY_TASK_ID} done  $(date)"
