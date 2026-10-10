# FEASIBILITY PROBE for the targeted large-n RATE run. Measures the two unknowns that decide
# whether the run is possible at all, and at what design.
#
# WHY. The headline ("not root-n over n in [300,5000]") is fitted from THREE sample sizes with a
# TWO-parameter log-log fit => ZERO residual df => it cannot distinguish a real non-root-n limit from
# a slow transient. Adding n = 20000 and 50000 gives 5 points / 3 residual df, enough to test
# CURVATURE. Everything else (mechanism hunting, the honest-limitation write-up) is downstream of
# that answer.
#
# THE RATE QUESTION DOES NOT NEED THE M4 BOOTSTRAP (77 min/rep at n=5000 -- the thing that makes
# the main run expensive). It needs only, per replicate:
#     point estimate MRD-hat(tau)   -> fold_sd across reps
#     M6 analytic cluster SE        -> mean_se across reps
# i.e. dat_sim + ste_prepare + .aalen_if + ste_estimate. No multiplier draws, no band, no bootstrap.
#
# UNKNOWN 1: seconds/rep at n = 5000 / 20000 / 50000 (aalen is NOT linear in n -- the event-time loop
#            is O(#event times x risk-set size)). This decides feasibility and whether ARC is needed.
# UNKNOWN 2: how many reps fold_sd needs. fold_sd = quantile(|est - mean|, 0.95)/qnorm(0.975); a 0.95
#            quantile is noisy, and we must resolve a CURVATURE, not just a slope.
#
#   Rscript dev/probe_largen_cost.R
suppressWarnings(suppressMessages({library(dplyr)}))
ROOT <- "."
for (f in list.files(file.path(ROOT, "steCI/R"), pattern = "[.]R$", full.names = TRUE)) source(f)

# One replicate, M6-only: exactly what the rate run needs and nothing more.
one_rep <- function(n, gamma.L, gamma.0, alpha.0, seed, t.hor = 1:5) {
  d  <- dat_sim(n = n, n.visit = 5, gamma.0 = gamma.0, gamma.L = gamma.L, alpha.0 = alpha.0,
                alpha.A = -0.04, alpha.L = 0.015, alpha.U = 0.015, seed = seed)
  ds <- ste_prepare(d, n.visit = 5, truncate = 0.95)
  fit <- .aalen_if(ds, t.hor = t.hor, weight_col = "ipw.stab.cum")
  # MRD(tau) by empirical standardization over the baseline-L population, + cluster-robust SE.
  L0 <- d$L[d$visit == 1]
  est <- se <- numeric(length(t.hor))
  for (s in seq_along(t.hor)) {
    B <- fit$B[s, ]
    s1 <- mean(exp(-(B[1] + B[2] + B[3] * L0)))
    s0 <- mean(exp(-(B[1] +        B[3] * L0)))
    est[s] <- s0 - s1
    g <- c(mean(exp(-(B[1] + B[3]*L0))) - mean(exp(-(B[1] + B[2] + B[3]*L0))),
           -mean(exp(-(B[1] + B[2] + B[3]*L0))) * -1,
           mean(L0 * exp(-(B[1] + B[3]*L0))) - mean(L0 * exp(-(B[1] + B[2] + B[3]*L0))))
    g <- c(s1 - s0, s1, mean(L0*exp(-(B[1]+B[2]+B[3]*L0))) - mean(L0*exp(-(B[1]+B[3]*L0))))
    ifc <- fit$IFc[, , s] %*% g
    se[s] <- sqrt(sum(ifc^2))
  }
  list(est = est, se = se, n_rows = nrow(ds))
}

cat("######## LARGE-n RATE RUN: FEASIBILITY PROBE ########\n\n")
cat("Per-rep work = dat_sim + ste_prepare + .aalen_if + standardization. NO bootstrap, NO band.\n\n")

cat("UNKNOWN 1: seconds per replicate vs n (the hard constraint)\n")
tm <- do.call(rbind, lapply(c(1000L, 5000L, 20000L), function(n) {
  t0 <- proc.time()[["elapsed"]]
  r <- one_rep(n, 3.0, -2.0, 0.25, 101)
  el <- proc.time()[["elapsed"]] - t0
  data.frame(n = n, secs_per_rep = round(el, 2), dat_seq_rows = r$n_rows,
             mrd_tau5 = round(r$est[5], 4), se_tau5 = round(r$se[5], 4))
}))
print(as.data.frame(tm))
sc <- lm(log(secs_per_rep) ~ log(n), data = tm)
pw <- unname(coef(sc)[2])
cat("\n  cost scaling: secs ~ n^", round(pw, 2), "\n", sep = "")
for (N in c(50000L)) {
  pred <- exp(predict(sc, newdata = data.frame(n = N)))
  cat("  EXTRAPOLATED n=", N, ": ", round(pred, 1), " s/rep  => 500 reps = ",
      round(pred * 500 / 3600, 1), " core-hours; on 32 cores = ",
      round(pred * 500 / 3600 / 32, 2), " h\n", sep = "")
}
cat("\n  BUDGET CHECK for the full design (5 n x 2 gamma.L x 2 alpha.0 x 2 gamma.0 x 500 reps):\n")
tot <- 0
for (N in c(300L, 1000L, 5000L, 20000L, 50000L)) {
  s <- if (N <= 20000) { i <- match(N, tm$n); if (!is.na(i)) tm$secs_per_rep[i] else exp(predict(sc, newdata = data.frame(n = N))) } else exp(predict(sc, newdata = data.frame(n = N)))
  tot <- tot + s * 500 * 8   # 8 = 2 gamma.L x 2 alpha.0 x 2 gamma.0
}
cat("    total ~", round(tot/3600, 1), "core-hours =", round(tot/3600/32, 1), "h on 32 cores\n")

cat("\nUNKNOWN 2: how many reps does fold_sd need? (MC error of the 0.95 quantile)\n")
cat("  Bootstrapping fold_sd from the EXISTING 1000-rep output (no new sims):\n")
rp <- readRDS(file.path(ROOT, "simulation/output/results/replicates_all.rds"))
pt <- rp %>% filter(method == "point", estimand == "mrd", tau == 5,
                    scenario == "n5000_evmed_cfstrong_pvlow", is.finite(estimate))
x <- pt$estimate
foldsd <- function(v) unname(quantile(abs(v - mean(v)), 0.95) / qnorm(0.975))
set.seed(1)
for (R in c(250, 500, 1000)) {
  bs <- replicate(400, foldsd(sample(x, R, replace = TRUE)))
  cat(sprintf("    R=%4d reps -> fold_sd = %.4f, MC SE = %.5f  (rel. MC SE = %.2f%%)\n",
              R, mean(bs), sd(bs), 100 * sd(bs) / mean(bs)))
}
cat("\n  Read: to resolve CURVATURE we must distinguish b=0.5 from b=0.25 across n.\n")
cat("  Over a 2.5x step in n, that is a log(fold_sd) difference of ~0.25*log(2.5) = 0.23 (~26%).\n")
cat("  A relative MC SE of ~3% per point is comfortably enough; ~10% would not be.\n")

cat("\nPROBE DONE\n")
