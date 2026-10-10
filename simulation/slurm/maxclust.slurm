#!/bin/bash
# =====================================================================
# maxclust.slurm -- max.clust sensitivity, UBC ARC (Sockeye).
# One array task = one (n=5000 cell x replicate-chunk) -> one RDS in simulation/output/maxclust/.
# Checkpointed: an existing output file is skipped, so a re-submit resumes in seconds.
#
# WHY. The published M3 cluster-robust SE uses timereg::aalen(..., clusters=id) at the DEFAULT
# max.clust = 1000, which quantile-bins the 5000 subject clusters at n=5000 (n=300/1000 are exact).
# The binned meat is UNBIASED for the exact subject-level meat, so the shift from pinning
# max.clust=NULL is noise about zero -- but that rests on a 20-replicate check (+0.37%, p=0.33).
# This tightens it: re-fit the cluster SE with and without the pin at every n=5000 cell, many reps.
#
# COMPUTES NO TRUTH. Emits per-replicate (gamma, se_default, se_pinned) via the SHIPPED cluster
# fit (byte-identical to ci_analytic.R); scoring is dev/maxclust_collect.R, at home.
#
#   SMOKE : bash simulation/slurm/submit_maxclust.sh smoke     # 4 tasks, minutes
#   FULL  : bash simulation/slurm/submit_maxclust.sh full      # 270 tasks
# =====================================================================
#SBATCH --account=YOUR_HPC_ACCOUNT
#SBATCH --job-name=steCI_maxclust
#SBATCH --array=1-540
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=0:45:00
# WALLTIME/ARRAY/MEM are set per-profile by submit_maxclust.sh, which overrides the values above.
#
# SIZING -- MEASURED on the ARC smoke (job 12346988): 29-47 s per aalen fit on the FAST evhigh cells,
#   peak RSS ~1.3 GB. FULL = 27 cells x 100 reps, 5 reps/chunk -> 540 tasks (10 fits/task); the evlow
#   cells carry ~1.8x the stacked rows and run slower/larger, so submit_maxclust.sh sets 0:45:00 / 6G.
#   The exact-subject (max.clust=NULL) fit is the slower of the two -- re-run the smoke if in doubt.
#
# MEMORY -- MEASURED: evhigh peaked ~1.3 GB on ARC; evlow (more rows) is larger, hence --mem=6G.
#SBATCH --output=simulation/output/logs/mc_%A_%a.out
#SBATCH --error=simulation/output/logs/mc_%A_%a.err
#SBATCH --mail-user=you@example.org
#SBATCH --mail-type=FAIL
# No --partition: the YOUR_HPC_ACCOUNT allocation routes to the right nodes.

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p simulation/output/maxclust simulation/output/logs

module load gcc/9.4.0
module load r/4.4.0
export R_LIBS="${R_LIBS:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"

export REPO_ROOT="$SLURM_SUBMIT_DIR"
export MC_OUT="$SLURM_SUBMIT_DIR/simulation/output/maxclust"
# MC_CELLS is COLON-separated (1-based indices into the 27 n=5000 cells; empty = all). sbatch
# --export is comma-separated, so a comma inside a value silently truncates the list.
export MC_CELLS="${MC_CELLS:-}"
export MC_REPS="${MC_REPS:-100}"
export MC_CHUNK="${MC_CHUNK:-10}"

echo "[mc] task ${SLURM_ARRAY_TASK_ID} start $(date) on $(hostname)"
Rscript dev/maxclust_arc.R
echo "[mc] task ${SLURM_ARRAY_TASK_ID} done  $(date)"
