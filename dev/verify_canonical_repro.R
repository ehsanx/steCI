# CLOSE THE REPRODUCTION QUESTION -- bit-for-bit, not statistically.
#
# The estimand panel's `adjustment` ablation measured +0.003792 +/- 0.000903 for
# n5000_evlow_cfweak_pvmed where the canonical run records +0.006896 +/- 0.000223 -- a 3.3 sigma
# gap it flagged as "a bigger problem than the ablation question". My earlier claim that the
# robustness re-run superseded this was WRONG: the robustness studies use their own self-contained
# fit_e1/fit2, not one_replicate.R -> ste_estimate. Different code path; it validates neither.
#
# This does not compare means. The canonical replicate seeds are exactly reproducible
# (run_sim.R:59: set.seed(SEED_BASE); all_seeds <- sample.int(.Machine$integer.max, N_REPS)), and
# replicates_all.rds holds every per-replicate point estimate. So we regenerate the seeds, rerun
# the PRODUCTION point-estimate path on the first few replicates, and compare estimate by
# estimate. A handful of exact matches settles it; a statistical comparison never could.
#
# The production path for the point estimate is one_replicate.R:28-53:
#     set.seed(rep_seed) -> dat_sim(...) -> ste_estimate(..., weight="stab")$gamma
# computed BEFORE any bootstrap, so B is irrelevant and this is cheap.
suppressMessages({library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang)})
for (f in list.files("steCI/R", pattern = "[.]R$", full.names = TRUE)) source(f)
source("simulation/scenarios.R")

SEED_BASE <- 20250228L; N_REPS <- 1000L
set.seed(SEED_BASE); all_seeds <- sample.int(.Machine$integer.max, N_REPS)

SC <- "n5000_evlow_cfweak_pvmed"
g  <- make_scenarios(); cell <- g[g$scenario == SC, ]
stopifnot(nrow(cell) == 1)
cat(sprintf("scenario %s: n=%d alpha.0=%s gamma.0=%s gamma.L=%s\n\n",
            SC, cell$n, format(cell$alpha.0), format(cell$gamma.0), format(cell$gamma.L)))

r <- readRDS("data/exact/results/replicates_all.rds")
can <- r[r$scenario == SC & r$method == "point" & r$estimand == "hd", c("rep", "estimate")]
can <- can[order(can$rep), ]
cat(sprintf("canonical: %d replicates, mean %.7f, sd %.7f\n",
            nrow(can), mean(can$estimate), sd(can$estimate)))

NREP <- as.integer(Sys.getenv("NREP", "12"))
cat(sprintf("\nreproducing the first %d replicates through the PRODUCTION path...\n\n", NREP))
cat(sprintf("%5s %14s %18s %18s %14s\n", "rep", "seed", "canonical", "reproduced", "difference"))
d <- numeric(NREP)
for (i in seq_len(NREP)) {
  set.seed(all_seeds[i])
  dd <- dat_sim(n = cell$n, n.visit = DGM_FIXED$n.visit,
                gamma.0 = cell$gamma.0, gamma.L = cell$gamma.L,
                alpha.0 = cell$alpha.0, alpha.A = DGM_FIXED$alpha.A,
                alpha.L = DGM_FIXED$alpha.L, alpha.U = DGM_FIXED$alpha.U)
  est <- ste_estimate(dd, n.visit = DGM_FIXED$n.visit, t.hor = T_HOR,
                      weight = "stab", estimand = "both")
  cv <- can$estimate[can$rep == i]
  d[i] <- est$gamma - cv
  cat(sprintf("%5d %14d %18.10f %18.10f %14.2e\n", i, all_seeds[i], cv, est$gamma, d[i]))
}
# Tolerance keyed to MAGNITUDE, not exact zero. Two different floating-point checks today
# demanded bit-identity and raised false alarms: identical arithmetic in a different evaluation
# order differs in the last ulp, so the floor is ~|x|*eps, not 0. Here |gamma| ~ 0.05 and
# eps = 2.2e-16, so the floor is ~1e-17 and anything under a few hundred ulp is "identical".
tol <- 1000 * max(abs(can$estimate[seq_len(NREP)])) * .Machine$double.eps
cat(sprintf("\nmax |difference| over %d replicates: %.3e   (tolerance %.3e = 1000 ulp)\n",
            NREP, max(abs(d)), tol))
cat(sprintf("  in units of one ulp on this magnitude: %.1f\n",
            max(abs(d)) / (max(abs(can$estimate[seq_len(NREP)])) * .Machine$double.eps)))
cat(if (max(abs(d)) < tol)
  "\n  => REPRODUCES EXACTLY, to last-ulp arithmetic. The canonical primary pipeline is sound.\n  The ablation agent's hand-rolled dat_sim -> ste_prepare -> aalen was the outlier: the\n  production path calls ste_estimate(), which is not the same object.\n" else
  "\n  => DOES NOT REPRODUCE. This is a real problem and outranks the ablation question.\n")
