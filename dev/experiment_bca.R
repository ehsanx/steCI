# Experiment: does a skewness+bias-corrected (BCa) bootstrap interval close the
# residual under-coverage of the marginal risk difference, and is the residual
# skewness (fixable) or bias (a point-estimator limitation)?
#
# Compares, per horizon, at n=300 (the STE-relevant small-n regime):
#   Wald (mrd +/- 1.96 SE), percentile, basic/pivot, and BCa.
# BCa uses z0 from the nonparametric bootstrap and acceleration a from the M6
# influence functions (a = (1/6) sum IF^3 / (sum IF^2)^{3/2}), avoiding a nested
# jackknife. Run in background:  Rscript dev/experiment_bca.R
suppressMessages({library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang)})
for (f in list.files("steCI/R", pattern = "[.]R$", full.names = TRUE)) source(f)
source("simulation/scenarios.R")
# TRUTHS: dev/truth_cache_exact.rds is the closed-form cache.
# This script previously read simulation/output/truth_cache.rds -- a SINGLE frozen
# Monte-Carlo draw, displaced from the exact values by +0.0012/+0.0023/+0.0031 -- so
# any coverage it printed was scored against the wrong target. Do not revert.
tv <- readRDS("dev/truth_cache_exact.rds")

TH <- 1:5; K <- 150L; B <- 200L; al <- 0.025
za <- stats::qnorm(al); z1a <- stats::qnorm(1 - al)
cells <- subset(make_scenarios(),
                n == 300 & prevalence == "med" & confounding == "moderate" &
                event_rate %in% c("low", "high"))
set.seed(99); seeds <- sample.int(1e8, K)

# Assert key-space coverage BEFORE the loop, so a missing truth fails at item 1
# rather than silently mis-scoring a cell.
{
  need <- unique(as.character(cells$alpha.0))
  miss <- setdiff(need, names(tv))
  if (length(miss)) stop("truth cache lacks alpha.0 key(s): ", paste(miss, collapse = ", "))
  stopifnot(all(vapply(tv[need], function(z) isTRUE(z$exact), logical(1))))
}

# point MRD, cluster SE, and M6 influence-function acceleration, per horizon
m6_bits <- function(d) {
  ds <- ste_prepare(d); fit <- .aalen_if(ds, TH); L0 <- d$L[d$visit == 1]
  nt <- length(TH); mrd <- se <- accel <- numeric(nt)
  for (s in seq_len(nt)) {
    B0 <- fit$B[s,1]; BA <- fit$B[s,2]; BL <- fit$B[s,3]
    S1 <- exp(-(B0+BA+BL*L0)); S0 <- exp(-(B0+BL*L0))
    mrd[s] <- mean(S0 - S1)
    g <- c(mean(S1-S0), mean(S1), mean(L0*(S1-S0)))
    ifc <- fit$IFc[,,s] %*% g
    se[s] <- sqrt(sum(ifc^2))
    accel[s] <- (sum(ifc^3)) / (6 * (sum(ifc^2))^1.5)
  }
  list(mrd = mrd, se = se, accel = accel)
}

for (ci in seq_len(nrow(cells))) {
  cell <- cells[ci,]; truth <- tv[[as.character(cell$alpha.0)]]$mrd
  meths <- c("wald","percentile","basic","bca")
  cov <- array(0, c(length(meths), length(TH))); nvalid <- 0L
  for (k in seq_len(K)) {
    d <- dat_sim(n = cell$n, gamma.0 = cell$gamma.0, gamma.L = cell$gamma.L,
                 alpha.0 = cell$alpha.0, seed = seeds[k])
    est <- tryCatch(m6_bits(d), error = function(e) NULL); if (is.null(est)) next
    resample <- .make_resampler(d)
    boot <- matrix(NA_real_, B, length(TH))
    for (b in seq_len(B)) {
      e <- tryCatch(ste_estimate(.resample_wrap <- resample(), estimand = "mrd", t.hor = TH),
                    error = function(e) NULL)
      if (!is.null(e)) boot[b,] <- e$mrd
    }
    nvalid <- nvalid + 1L
    for (s in seq_len(length(TH))) {
      bt <- boot[is.finite(boot[,s]), s]; if (length(bt) < 50) next
      pt <- est$mrd[s]; se <- est$se[s]; a <- est$accel[s]
      z0 <- stats::qnorm(pmin(pmax(mean(bt < pt), 1/(2*B)), 1 - 1/(2*B)))
      a1 <- stats::pnorm(z0 + (z0 + za)/(1 - a*(z0 + za)))
      a2 <- stats::pnorm(z0 + (z0 + z1a)/(1 - a*(z0 + z1a)))
      q  <- stats::quantile(bt, c(al, 1-al, a1, a2), names = FALSE)
      lo <- c(wald = pt - z1a*se, percentile = q[1], basic = 2*pt - q[2], bca = q[3])
      hi <- c(wald = pt + z1a*se, percentile = q[2], basic = 2*pt - q[1], bca = q[4])
      cov[, s] <- cov[, s] + as.numeric(lo <= truth[s] & hi >= truth[s])
    }
  }
  cov <- cov / nvalid
  cat(sprintf("\n=== cell %s (n=300, event=%s), reps=%d, B=%d ===\n",
              cell$scenario, cell$event_rate, nvalid, B))
  dimnames(cov) <- list(meths, paste0("t", TH))
  print(round(cov, 3))
  cat("mean over horizons:", paste(sprintf("%s=%.3f", meths, rowMeans(cov)), collapse="  "), "\n")
  cat("late (t4-5):       ", paste(sprintf("%s=%.3f", meths, rowMeans(cov[,4:5])), collapse="  "), "\n")
}
cat("\nEXPERIMENT_BCA DONE\n")
