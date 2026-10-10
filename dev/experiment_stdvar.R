# Focused E2 coverage experiment for the random-standardisation (unconditional)
# MRD variance. Computes the weighted-Aalen influence functions ONCE per replicate
# (via .aalen_if) and forms, at each horizon, the analytic cluster-robust Wald
# interval and the multiplier interval under BOTH standardisation="fixed"
# (conditional; the manuscript's primary M6) and "random" (unconditional; adds the
# g-computation empirical-process term). Coverage is against the cached true MRD.
#
# The nonparametric bootstrap (M4) is NOT re-run per replicate here (it refits the
# whole pipeline and is infeasible off-cluster); dev/validate_stdvar.R shows that
# the random-standardisation SE matches the bootstrap SE on a single dataset.
#   Rscript dev/experiment_stdvar.R
suppressMessages({library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang)})
for (f in list.files("steCI/R", pattern = "[.]R$", full.names = TRUE)) source(f)
# TRUTHS: dev/truth_cache_exact.rds is the closed-form cache (dev/truth_closed_form.R).
# An earlier version read simulation/output/truth_cache.rds -- a SINGLE frozen Monte-Carlo
# draw that simulation/metrics.R flags as SUPERSEDED (its E2/MRD targets are displaced by up
# to +0.0016) -- so any coverage it printed was scored against the wrong target. Same fix as
# dev/experiment_bca.R. Do not revert.
tv <- readRDS("dev/truth_cache_exact.rds")
TH <- 1:5; z <- stats::qnorm(0.975); a <- 0.025

# representative grid (medium prevalence/confounding); n=5000 at medium event only
cells <- data.frame(
  n       = c(300L, 300L, 300L, 1000L, 1000L, 1000L, 5000L),
  alpha.0 = c(0.15, 0.25, 0.40, 0.15,  0.25,  0.40,  0.25),
  reps    = c(400L, 400L, 400L, 400L,  400L,  400L,  200L),
  Bmult   = 200L, gamma.0 = -1, gamma.L = 1.5)
cells$ev <- c("low", "med", "high")[match(cells$alpha.0, c(0.15, 0.25, 0.40))]

# one-replicate M6: fit once, return cluster-Wald & multiplier coverage indicators
one <- function(d, truth, Bmult, seed) {
  ds <- ste_prepare(d); fit <- .aalen_if(ds, TH)
  L0 <- d$L[d$visit == 1]
  Lc <- ds$L.baseline[match(fit$clusters, as.character(ds$id))]
  K  <- dim(fit$IFc)[1]; nt <- length(TH)
  set.seed(seed); Xi <- matrix(sample(c(-1, 1), K * Bmult, replace = TRUE), K, Bmult)
  meth <- c("cl_fixed", "cl_random", "mult_fixed", "mult_random")
  hit <- matrix(0, length(meth), nt, dimnames = list(meth, NULL))
  for (s in seq_len(nt)) {
    B0 <- fit$B[s, 1]; BA <- fit$B[s, 2]; BL <- fit$B[s, 3]
    S1 <- exp(-(B0 + BA + BL * L0)); S0 <- exp(-(B0 + BL * L0))
    mrd <- mean(S0 - S1)
    g   <- c(mean(S1 - S0), mean(S1), mean(L0 * (S1 - S0)))
    base <- as.vector(fit$IFc[, , s] %*% g)                 # fixed IF
    S1c <- exp(-(B0 + BA + BL * Lc)); S0c <- exp(-(B0 + BL * Lc))
    stdv <- ((S0c - S1c) - mean(S0c - S1c)) / K             # g-computation term (1/K)
    rand <- base + stdv                                     # random IF
    se_f <- sqrt(sum(base^2)); se_r <- sqrt(sum(rand^2))
    # multiplier draws
    mf <- mrd + as.vector(crossprod(base, Xi))
    mr <- mrd + as.vector(crossprod(rand, Xi))
    qf <- stats::quantile(mf, c(a, 1 - a), names = FALSE)
    qr <- stats::quantile(mr, c(a, 1 - a), names = FALSE)
    tr <- truth[s]
    hit["cl_fixed",    s] <- (mrd - z*se_f <= tr && tr <= mrd + z*se_f)
    hit["cl_random",   s] <- (mrd - z*se_r <= tr && tr <= mrd + z*se_r)
    hit["mult_fixed",  s] <- (qf[1] <= tr && tr <= qf[2])
    hit["mult_random", s] <- (qr[1] <= tr && tr <= qr[2])
  }
  hit
}

summ <- list()
for (ci in seq_len(nrow(cells))) {
  cell <- cells[ci, ]; truth <- tv[[as.character(cell$alpha.0)]]$mrd
  meth <- c("cl_fixed", "cl_random", "mult_fixed", "mult_random")
  cov <- matrix(0, length(meth), length(TH), dimnames = list(meth, NULL))
  nv <- 0L
  for (r in seq_len(cell$reps)) {
    seed <- 1000L * ci + r
    d <- tryCatch(dat_sim(n = cell$n, alpha.0 = cell$alpha.0, gamma.0 = cell$gamma.0,
                          gamma.L = cell$gamma.L, seed = seed), error = function(e) NULL)
    if (is.null(d)) next
    h <- tryCatch(one(d, truth, cell$Bmult, seed), error = function(e) NULL)
    if (is.null(h)) next
    nv <- nv + 1L; cov <- cov + h
  }
  cov <- cov / nv
  cat(sprintf("\n=== n=%d ev=%s (alpha0=%.2f) reps=%d ===\n", cell$n, cell$ev, cell$alpha.0, nv))
  print(round(cov, 3))
  cat("mean/horizon:", paste(sprintf("%s=%.3f", meth, rowMeans(cov)), collapse = "  "), "\n")
  summ[[ci]] <- data.frame(n = cell$n, ev = cell$ev, t(round(rowMeans(cov), 3)))
}
cat("\n=== SUMMARY (mean coverage over horizons) ===\n")
print(do.call(rbind, summ))
cat("\nEXPERIMENT_STDVAR DONE\n")
