# STEP 1: reconcile 0.678 (localisation run) vs 0.726 (ARC truncation run) -- same cell, same
# quantity (strong_conf, n=5000, tau=5, Jensen-free rms_se/emp_sd), ~7% apart.
#
# TWO QUESTIONS, both answered by measurement rather than argument:
#
# (A) DO THE TWO CODE PATHS AGREE? Reading the source they should be identical:
#     ci_lef.R:84   g <- c(mean(S1 - S0), mean(S1), mean(L0 * (S1 - S0)))
#     localisation  g <- c(mean(S1)-mean(S0), mean(S1), mean(L0*S1)-mean(L0*S0))   <- expands to the same
#     Both use L0 = dat.long$L[visit==1], ste_prepare(truncate=0.95), weight_col="ipw.stab.cum",
#     and se_cl = sqrt(sum((IFc %*% g)^2)). So on the SAME SEED they must agree to machine precision.
#     If they do not, one path has a bug and everything built on it is suspect.
#
# (B) IS R=300 EVEN ENOUGH FOR THE LOCALISATION TEST? The plan review measured emp_sd's relative MC
#     error at 17-19% at R=500 in these high-kurtosis cells (kurtosis ~60 at alpha.0=0.25/gamma.L=3/
#     tau=5/n=5000), which implies ~22% at R=300. If so the localisation LEVELS (0.682/0.689/0.681)
#     are far too noisy to quote -- BUT the per-component COMPARISON is paired (same replicates), so
#     the DIFFERENCES may still be precise, and the claim was about differences (B_L 0.986 vs B0 0.804).
#     Bootstrap both, and report the levels and the differences separately.
#
#   Rscript dev/reconcile_paths.R
suppressWarnings(suppressMessages({library(dplyr); library(parallel)}))
ROOT <- "."
NREC <- as.integer(Sys.getenv("REC_N", "40"))    # seeds for the path comparison (exactness, not means)
TAU <- 1:5

## ---------------------------------------------------------------- (A) path agreement
cat("######## (A) DO THE TWO PATHS AGREE ON THE SAME SEED? ########\n\n")
for (f in list.files(file.path(ROOT, "steCI/R"), pattern = "[.]R$", full.names = TRUE)) source(f)

path_cilef <- function(d, seed) {
  lf <- ci_lef(d, t.hor = TAU, weight = "stab", B = 50, approach = "a1", seed = seed)
  m <- lf$mrd[lf$mrd$type == "analytic_cluster", ]
  data.frame(tau = m$tau, est = m$estimate, se = m$se)
}
path_hand <- function(d, seed) {
  ds <- ste_prepare(d, n.visit = 5, truncate = 0.95)
  fit <- .aalen_if(ds, t.hor = TAU, weight_col = "ipw.stab.cum")
  L0 <- d$L[d$visit == 1]
  out <- lapply(seq_along(TAU), function(s) {
    B <- fit$B[s, ]
    S1 <- exp(-(B[1] + B[2] + B[3] * L0)); S0 <- exp(-(B[1] + B[3] * L0))
    g <- c(mean(S1) - mean(S0), mean(S1), mean(L0 * S1) - mean(L0 * S0))
    ifm <- fit$IFc[, , s] %*% g
    data.frame(tau = TAU[s], est = mean(S0) - mean(S1), se = sqrt(sum(ifm^2)))
  })
  do.call(rbind, out)
}

cmp <- do.call(rbind, lapply(1:NREC, function(i) {
  sd_ <- 9e6 + i
  d <- dat_sim(n = 1000L, n.visit = 5, gamma.0 = -2, gamma.L = 3, alpha.0 = 0.25,
               alpha.A = -0.04, alpha.L = 0.015, alpha.U = 0.015, seed = sd_)
  a <- path_cilef(d, sd_); b <- path_hand(d, sd_)
  data.frame(seed = sd_, tau = a$tau, d_est = a$est - b$est, d_se = a$se - b$se)
}))
cat("  seeds:", NREC, "| max |est_cilef - est_hand| =", format(max(abs(cmp$d_est)), digits = 3), "\n")
cat("                  max |se_cilef  - se_hand |  =", format(max(abs(cmp$d_se)), digits = 3), "\n")
cat("  VERDICT:", if (max(abs(cmp$d_est)) < 1e-10 && max(abs(cmp$d_se)) < 1e-10)
      "PATHS ARE IDENTICAL -> the 0.678 vs 0.726 gap is Monte Carlo, not a bug."
    else "*** PATHS DIFFER -- there is a bug. Investigate before trusting either. ***", "\n")

## ---------------------------------------------------------------- (B) is R=300 enough?
cat("\n######## (B) IS R=300 ENOUGH FOR THE LOCALISATION TEST? ########\n\n")
r <- readRDS(file.path(ROOT, "dev/localisation_results.rds"))
w <- r %>% filter(cell == "strong_conf", n == 5000, tau == 5)
cat("  cell: strong_conf, n=5000, tau=5 |", n_distinct(w$seed), "reps\n")
cat("  excess kurtosis of the point estimates by component:\n")
for (q in c("B0", "B_A", "B_L", "MRD")) {
  x <- w$est[w$quantity == q]; z <- x - mean(x)
  cat(sprintf("    %-4s kurt = %6.1f\n", q, mean(z^4)/mean(z^2)^2 - 3))
}

ratio_rms <- function(df) sqrt(mean(df$se^2)) / sd(df$est)
set.seed(1)
sds <- unique(w$seed); NB <- 600
bs <- replicate(NB, {
  ix <- sample(sds, length(sds), replace = TRUE)
  vapply(c("B0", "B_A", "B_L", "MRD"), function(q) {
    sub <- w[w$quantity == q, ]
    ratio_rms(sub[match(ix, sub$seed), ])
  }, numeric(1))
})
cat("\n  BOOTSTRAPPED LEVELS (600 resamples over replicates):\n")
for (q in c("B0", "B_A", "B_L", "MRD")) {
  v <- bs[q, ]
  cat(sprintf("    %-4s ratio = %.3f  MC SE = %.3f  (rel %.1f%%)  95%% CI [%.3f, %.3f]\n",
              q, mean(v), sd(v), 100*sd(v)/mean(v), quantile(v, .025), quantile(v, .975)))
}
cat("\n  *** THE QUESTION THAT MATTERS *** -- the claim was about DIFFERENCES (B_L 0.986 vs B0 0.804).\n")
cat("  Differences are PAIRED (same replicates), so they can be precise even when levels are not:\n")
for (p in list(c("B_L","B0"), c("B_L","B_A"), c("B_A","B0"))) {
  v <- bs[p[1], ] - bs[p[2], ]
  cat(sprintf("    %-4s - %-4s = %+.4f  MC SE = %.4f  z = %+5.2f   %s\n",
              p[1], p[2], mean(v), sd(v), mean(v)/sd(v),
              if (abs(mean(v)/sd(v)) > 1.96) "DIFFERENT" else "NOT DISTINGUISHABLE"))
}
cat("\n  Claimed difference B_L - B0 = 0.986 - 0.804 = +0.182. Is that excluded?\n")
v <- bs["B_L", ] - bs["B0", ]
cat(sprintf("    observed %+.4f, 95%% CI [%+.4f, %+.4f] -> the claimed +0.182 is %s\n",
            mean(v), quantile(v, .025), quantile(v, .975),
            if (0.182 < quantile(v, .025) || 0.182 > quantile(v, .975)) "EXCLUDED" else "NOT excluded"))
cat("\nRECONCILE DONE\n")
