#!/bin/bash
# =====================================================================
# aggregate.slurm -- combine all per-task RDS into metric tables.
# Submit AFTER sim.slurm (see submit_all.sh for the dependency wiring).
# =====================================================================
#SBATCH --account=YOUR_HPC_ACCOUNT
#SBATCH --job-name=steCI_agg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=simulation/output/logs/agg_%j.out
#SBATCH --error=simulation/output/logs/agg_%j.err
#SBATCH --mail-user=you@example.org
#SBATCH --mail-type=END,FAIL

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
module load gcc/9.4.0
module load r/4.4.0
export R_LIBS="${R_LIBS:-$HOME/R/x86_64-pc-linux-gnu-library/4.4}"
export REPO_ROOT="$SLURM_SUBMIT_DIR"
export SIM_OUT="$SLURM_SUBMIT_DIR/simulation/output/intermediate"
export SIM_RESULTS="$SLURM_SUBMIT_DIR/simulation/output/results"
export SIM_TRUTH_CACHE="${SIM_TRUTH_CACHE:-$SLURM_SUBMIT_DIR/dev/truth_cache_exact.rds}"
export SIM_PSEUDO_N="${SIM_PSEUDO_N:-500000}"

Rscript simulation/aggregate.R
