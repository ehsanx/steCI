#!/bin/bash
# =====================================================================
# sim.slurm -- steCI ADEMP simulation, UBC ARC (Sockeye) SLURM ARRAY.
# One array task = one (scenario x rep-chunk) cell -> one RDS.
# Copy the repo to your cluster scratch space (e.g. $SCRATCH/steCI/) and
# submit from there. NEVER run Rscript on the login nodes.
#
#   Use submit_all.sh -- it sets the array size and CHUNK together and cannot get them out of step.
#   SMOKE : bash simulation/slurm/submit_all.sh smoke
#   FULL  : bash simulation/slurm/submit_all.sh full
#
# Array size = (n_scenarios = 81) x ceil(SIM_N_REPS / SIM_CHUNK).
#   FULL: 81 x ceil(1000/32) = 81 x 32 = 2592  -> --array=1-2592
# =====================================================================
#SBATCH --account=YOUR_HPC_ACCOUNT
#SBATCH --job-name=steCI_sim
#SBATCH --array=1-2592
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
# SIM_CHUNK MUST BE <= --cpus-per-task (32). This is not advice: CHUNK=100 on 32 cores is exactly what
# cost the 2026-07 run 10 of 81 cells at walltime -- three landed on 640 reps = 2x32 mclapply waves
# killed partway through the third -- against this file's own earlier warning. submit_all.sh now sets
# CHUNK=32 = ONE wave per task. n=5000 is the heavy case: the M4 nonparametric bootstrap refits the
# whole pipeline B times (~77 min/rep), so one wave of 32 reps on 32 cores is ~77 min -- inside 12 h
# with a ~9x margin, whereas 100 reps needs 4 waves and blows it.
#SBATCH --output=simulation/output/logs/sim_%A_%a.out
#SBATCH --error=simulation/output/logs/sim_%A_%a.err
#SBATCH --mail-user=you@example.org
#SBATCH --mail-type=FAIL
# No --partition: the YOUR_HPC_ACCOUNT allocation routes to the right nodes.

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p simulation/output/intermediate simulation/output/logs

module load gcc/9.4.0
module load r/4.4.0
# VERIFY on first login: `echo $HOME` then `ls $HOME/R/x86_64-pc-linux-gnu-library/4.4`.
# Install deps once interactively: install.packages(c("dplyr","tidyr","timereg","survival","rlang"))
export R_LIBS="${R_LIBS:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"

export REPO_ROOT="$SLURM_SUBMIT_DIR"
export SIM_OUT="$SLURM_SUBMIT_DIR/simulation/output/intermediate"
export SIM_N_REPS="${SIM_N_REPS:-1000}"
export SIM_CHUNK="${SIM_CHUNK:-32}"
export SIM_B="${SIM_B:-200}"
export SIM_SEED_BASE="${SIM_SEED_BASE:-20250228}"
export SIM_JACK_MAXN="${SIM_JACK_MAXN:-0}"   # 0 = M5 jackknife off, as in submit_all.sh and the canonical run

echo "[sim] task ${SLURM_ARRAY_TASK_ID} start $(date) on $(hostname) cores=${SLURM_CPUS_PER_TASK}"
Rscript simulation/run_sim.R
echo "[sim] task ${SLURM_ARRAY_TASK_ID} done  $(date)"
