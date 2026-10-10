# Validation of the M6 machinery (.aalen_if + ci_lef) against ground truth.
# (a) the custom weighted Aalen fit reproduces timereg::aalen cumulative
#     coefficients to ~machine precision; (b) the M6 cluster-robust analytic and
#     multiplier SEs for the marginal risk difference agree with the
#     nonparametric bootstrap SE. Backs the "M6 matches the bootstrap" claim.
#   Rscript dev/validate_m6.R
suppressMessages({library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang)})
for (f in list.files("steCI/R", pattern = "[.]R$", full.names = TRUE)) source(f)

cat("== (a) .aalen_if vs timereg::aalen cumulative coefficients ==\n")
set.seed(101); d <- dat_sim(n = 800, seed = 101); ds <- ste_prepare(d)
th <- c(1, 2, 3, 4)
fit <- .aalen_if(ds, t.hor = th, weight_col = "ipw.stab.cum")
ref <- aalen(Surv(time.new, time.stop.new, event) ~ A.baseline + L.baseline,
             data = ds, n.sim = 0, weights = ds$ipw.stab.cum)
cum <- ref$cum
for (s in seq_along(th)) {
  ic <- max(which(cum[, 1] <= th[s]))
  cat(sprintf("  t=%d  max|B_custom - B_timereg| = %.1e\n", th[s],
              max(abs(fit$B[s, ] - cum[ic, 2:4]))))
}

cat("\n== (b) M6 SE vs nonparametric bootstrap SE for the MRD (n=1000) ==\n")
set.seed(202); d2 <- dat_sim(n = 1000, seed = 202)
lf <- ci_lef(d2, t.hor = 1:5, B = 1000, seed = 5)
bo <- ci_bootstrap(d2, t.hor = 1:5, estimand = "mrd", B = 300, seed = 7)
bse <- tapply(bo$mrd$se, bo$mrd$tau, `[`, 1)
acl <- lf$mrd$se[lf$mrd$type == "analytic_cluster"]
mul <- lf$mrd$se[lf$mrd$type == "multiplier"]
arow <- lf$mrd$se[lf$mrd$type == "analytic_row"]
cat("  t : bootstrap | M6 analytic_cluster | M6 multiplier | M6 analytic_row\n")
for (k in 1:5)
  cat(sprintf("  %d : %.4f    | %.4f            | %.4f       | %.4f\n",
              k, bse[k], acl[k], mul[k], arow[k]))
cat("\nExpected: cluster/multiplier track the bootstrap; row-level is too small at late t.\n")
