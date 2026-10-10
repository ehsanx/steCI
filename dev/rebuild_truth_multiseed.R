# The single-frozen-draw truth defect: the E1 target `gamma_true` was ONE frozen Monte-Carlo draw, whose replication SD
# is comparable to the E1 bias it is used to measure. Being a single draw shared across every
# replicate and scenario it is a FIXED displacement that does not average out, and mc_se_cov
# models only binomial replicate noise.
#
# ROUTE. The first suggested fix -- raise n_hd from 3e4 to 5e5 -- is INFEASIBLE and
# also inefficient, and both are now measured rather than assumed:
#   cost of one fit at n_hd=3e4 is 88s / 130s / 165s at alpha.0 = 0.15 / 0.25 / 0.4, and
#   timereg::aalen() is ~quadratic in n (a 3x3 solve per event time over the risk set).
#   Averaging K seeds:  var ∝ 1/K,     cost ∝ K      -> var x cost is CONSTANT.
#   Raising n_hd:       var ∝ 1/n_hd,  cost ∝ n_hd^2 -> var x cost ∝ n_hd, strictly WORSE.
# So we average independent draws at the current n_hd. This also yields the reportable MC SE for
# the target, which raising n_hd would not.
#
# WHY EACH DRAW IS INDEPENDENT (and why none of them reproduces the cached value). true_values()
# does U <- rnorm(n); L0 <- rnorm(n, U, 1), so each vector draw consumes n numbers and the RNG
# stream DESYNCHRONISES with n: for a common seed, U agrees on the overlap but L0 does not
# (verified). A run at n = 3e4 is therefore an independent realisation, not a sub-sample of the
# canonical n = 5e5 run. That is exactly what we want for averaging -- but it means the cached
# value is one more draw from this same distribution, not a reference point any run should match.
#
# E2 IS UNAFFECTED. The survival truths s1/s0/mrd come from KM on the FULL n_pseudo = 5e5 under
# the canonical seed; gamma_true is fitted on only the first n_hd rows (true_values.R:60). This
# script touches gamma_true ONLY; dev/apply_m8_truth.R asserts the E2 entries are byte-identical.
#
#   Rscript dev/rebuild_truth_multiseed.R          # M8_NSEED (default 150), M8_CORES (default 11)
suppressMessages({library(parallel)})

NSEED <- as.integer(Sys.getenv("M8_NSEED", "150"))
CORES <- as.integer(Sys.getenv("M8_CORES", "11"))
N_HD  <- 3e4L
BASE  <- 77000000L   # NOT the canonical 20222022: see the note above -- these are fresh draws.
CANON <- "data/truth_cache.rds"
ROOT  <- normalizePath(".", winslash = "/")

stopifnot(file.exists(CANON))
tc  <- readRDS(CANON)
a0s <- names(tc)

jobs <- expand.grid(a0 = a0s, k = seq_len(NSEED), stringsAsFactors = FALSE)
jobs$seed <- BASE + seq_len(nrow(jobs))          # a distinct seed per (level, draw)
cat(sprintf("multi-seed truth rebuild\n  levels: %s\n  seeds/level: %d  jobs: %d  cores: %d  n_hd: %d\n",
            paste(a0s, collapse = ", "), NSEED, nrow(jobs), CORES, N_HD))
cat(sprintf("  est. serial CPU: %.1f h  -> est. wall: %.0f min\n\n",
            NSEED * 383 / 3600, NSEED * 383 / 60 / CORES))

cl <- makePSOCKcluster(CORES)
on.exit(stopCluster(cl), add = TRUE)
# PSOCK workers inherit NOTHING -- every object and every source() must be set up worker-side.
clusterExport(cl, c("ROOT", "N_HD", "jobs"), envir = environment())
invisible(clusterEvalQ(cl, {
  setwd(ROOT)
  suppressMessages({library(timereg); library(survival)})
  for (f in list.files("steCI/R", pattern = "[.]R$", full.names = TRUE)) source(f)
  source("simulation/scenarios.R")
  TRUE
}))

t0 <- Sys.time()
out <- parLapply(cl, seq_len(nrow(jobs)), function(i) {
  a0 <- jobs$a0[i]; sd_ <- jobs$seed[i]
  tv <- true_values(n = N_HD, n.visit = DGM_FIXED$n.visit, t.hor = T_HOR,
                    alpha.0 = as.numeric(a0), alpha.A = DGM_FIXED$alpha.A,
                    alpha.L = DGM_FIXED$alpha.L, alpha.U = DGM_FIXED$alpha.U,
                    seed = sd_, n_hd = N_HD)
  list(a0 = a0, seed = sd_, gamma = tv$gamma_true)
}, chunk.size = 1)
cat(sprintf("done in %.1f min\n\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

res <- list()
for (a0 in a0s) {
  g <- vapply(Filter(function(z) z$a0 == a0, out), function(z) z$gamma, numeric(1))
  g <- g[is.finite(g)]
  res[[a0]] <- list(gamma = g, n = length(g), mean = mean(g), sd = sd(g),
                    mcse = sd(g) / sqrt(length(g)), cached = tc[[a0]]$gamma_true)
}

cat("================ RESULT ================\n")
cat(sprintf("%-7s %4s %11s %11s %11s %11s %11s\n",
            "alpha.0", "K", "cached", "mean", "sd(1 draw)", "MCSE(mean)", "shift"))
for (a0 in a0s) with(res[[a0]],
  cat(sprintf("%-7s %4d %11.6f %11.6f %11.6f %11.6f %+11.6f\n",
              a0, n, cached, mean, sd, mcse, mean - cached)))
cat("\nReported E1 |bias| is ~0.005-0.0065, so compare MCSE(mean) against that.\n")
cat("'shift' is cached-minus-new: how far the frozen single draw sat from the target.\n")
cat(sprintf("shift/MCSE (how many SEs the old draw was off): %s\n",
            paste(sprintf("%s=%.1f", a0s,
              vapply(a0s, function(a) abs(res[[a]]$mean - res[[a]]$cached) / res[[a]]$mcse, 0)),
              collapse = "  ")))

saveRDS(res, "dev/truth_multiseed.rds")
cat("\nwrote dev/truth_multiseed.rds\n")
