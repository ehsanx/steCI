#!/usr/bin/env bash
# Prove the tie switch / time.grid leave the continuous-time simulation untouched (package
# defaults): (1) three complete canonical replicates, every method, vs replicates_all.rds;
# (2) dev/verify_canonical_repro.R, the 12 first n = 5000 point estimates.
# Optional: HEAD_R_DIR=<a copy of the pre-change steCI/R> also runs replicate (1a) on it.
# Run from the repository root:  bash dev/heart_ties/run_canonical_checks.sh
set -u
R=dev/heart_ties/runs
mkdir -p "$R"
Rscript dev/heart_ties/verify_default_path.R n300_evlow_cfweak_pvhigh 1   > "$R/verify_default_path_1.log" 2>&1 &
Rscript dev/heart_ties/verify_default_path.R n300_evhigh_cfstrong_pvlow 2 > "$R/verify_default_path_2.log" 2>&1 &
Rscript dev/heart_ties/verify_default_path.R n1000_evmed_cfmoderate_pvmed 3 > "$R/verify_default_path_3.log" 2>&1 &
NREP=12 Rscript dev/verify_canonical_repro.R > "$R/verify_canonical_repro.log" 2>&1 &
if [ -n "${HEAD_R_DIR:-}" ]; then
  STECI_R_DIR="$HEAD_R_DIR" Rscript dev/heart_ties/verify_default_path.R n300_evlow_cfweak_pvhigh 1 \
    > "$R/verify_default_path_1_HEADcode.log" 2>&1 &
fi
wait
echo CANONICAL_CHECKS_DONE
