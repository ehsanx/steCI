#!/bin/bash
# =====================================================================
# robust.slurm -- re-run the three appendix robustness studies, UBC ARC (Sockeye).
# One array task = one (study x cell x replicate-chunk) -> one RDS in
# simulation/output/robust/. Checkpointed: an existing output file is skipped, so a
# re-submit resumes in seconds.
#
# WHY THIS EXISTS. All three studies score raw coverage against ONE frozen Monte-Carlo
# draw of their own truth (the single-frozen-draw truth defect, applied to the robustness appendix). The truths
# are now exact for all three (dev/robustness_truth_exact.rds, dev/truth_misspec_exact.rds,
# dev/truth_tvfx_exact.rds), but the studies SAVED ONLY 9 SUMMARY ROWS EACH and
# discarded the per-replicate estimates -- so they cannot be re-scored, only re-run.
#
# THIS JOB COMPUTES NO TRUTH. It emits per-replicate (gamma_hat, se_model, se_row, se_cl)
# and nothing else; scoring happens on the workstation in dev/robust_collect.R, so any
# truth revision is a seconds-long re-score rather than another cluster job.
#
# Seeding is preserved from the original scripts (seed = 1000*gi + r), so this reproduces
# the SAME replicates and every coverage change is attributable to the truth alone.
#
#   SMOKE : bash simulation/slurm/submit_robust.sh smoke     # 9 tasks, minutes
#   FULL  : bash simulation/slurm/submit_robust.sh full      # 357 tasks
# =====================================================================
#SBATCH --account=YOUR_HPC_ACCOUNT
#SBATCH --job-name=steCI_robust
#SBATCH --array=1-357
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=0:45:00
# WALLTIME/ARRAY are set per-profile by submit_robust.sh, which overrides the values above.
#
# SIZING -- now MEASURED ON ARC (smoke, 24 tasks, 2026-07-19), not extrapolated.
# Per-replicate wall cost on ARC: 0.9 s at n=300, 5.0 s at n=1000, 64.5 s at n=5000.
# ARC is ~2.1x slower per core than the author's workstation (0.84 / 2.50 / 30.4 s), NOT the ~4x
# the truth run saw -- that factor came from a different workload (aalen on a 6e4-row pseudo-RCT
# rather than on stacked sequential trials), so the original budget was conservative.
# With chunk sizes 250 / 100 / 10 replicates by n (dev/robust_arc.R), that gives per-task walltime:
#     n=300   3.8 min      n=1000   8.3 min      n=5000  10.8 min
# Total compute across all 357 tasks: ~60 core-hours.
# 0:45:00 is ~4x headroom on the worst task. Kept well above the measured 10.8 min because these
# are single-core tasks that Slurm may pack densely, and the truth run's lesson was that
# concurrent fits contend for memory bandwidth. Do NOT inflate further: Slurm backfills SHORT
# jobs into gaps, so an over-long request queues longer than it runs.
#
# MEMORY: peak RSS MEASURED ON ARC at 455 MB (worst task, n=5000). The earlier 8G was sized off a
# 229 MB local R-heap figure plus guesswork about what timereg/BLAS allocate outside the heap;
# the measurement says 2G is already 4.4x headroom. Smaller tasks backfill faster, which is the
# whole reason this is 357 small tasks rather than a few fat ones.
#
# Do NOT inflate --time blindly: Slurm backfills SHORT jobs into gaps, so an over-long request
# queues far longer than it needs to.
#SBATCH --output=simulation/output/logs/robust_%A_%a.out
#SBATCH --error=simulation/output/logs/robust_%A_%a.err
#SBATCH --mail-user=you@example.org
#SBATCH --mail-type=FAIL
# No --partition: the YOUR_HPC_ACCOUNT allocation routes to the right nodes.

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p simulation/output/robust simulation/output/logs

module load gcc/9.4.0
module load r/4.4.0
export R_LIBS="${R_LIBS:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"

export REPO_ROOT="$SLURM_SUBMIT_DIR"
export ROBUST_OUT="$SLURM_SUBMIT_DIR/simulation/output/robust"
# ROBUST_STUDIES is COLON-separated. sbatch --export is itself comma-separated, so a comma inside
# a value silently truncates the list -- that bug cost the truth run its first smoke submission,
# where the n_hd sweep never ran and 6 of 9 tasks were no-ops. Commas still work for local use.
export ROBUST_STUDIES="${ROBUST_STUDIES:-2confounder:misspec:tvfx}"

echo "[robust] task ${SLURM_ARRAY_TASK_ID} start $(date) on $(hostname)"
Rscript dev/robust_arc.R
echo "[robust] task ${SLURM_ARRAY_TASK_ID} done  $(date)"
