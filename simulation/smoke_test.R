# Smoke test: exercise the steCI pipeline end-to-end on tiny inputs.
# Sources the package R files directly (no install needed).
# Run:  Rscript simulation/smoke_test.R
suppressMessages({
  library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang)
})

pkg_dir <- file.path("steCI", "R")
if (!dir.exists(pkg_dir)) pkg_dir <- file.path("..", "steCI", "R")
for (f in list.files(pkg_dir, pattern = "\\.R$", full.names = TRUE)) source(f)

cat("==== dat_sim ====\n")
d <- dat_sim(n = 300, seed = 1)
cat("rows:", nrow(d), " subjects:", length(unique(d$id)),
    " events:", sum(d$event), " neg_haz_frac:", round(attr(d, "neg_haz_frac"), 4), "\n")

cat("\n==== ste_estimate (both estimands) ====\n")
est <- ste_estimate(d, estimand = "both", t.hor = 1:5)
cat("gamma (estimated marginal constant hazard difference):", round(est$gamma, 4), "\n")
print(round(rbind(surv1 = est$surv1, surv0 = est$surv0,
                  survdiff = est$survdiff, mrd = est$mrd), 4))

cat("\n==== ci_analytic (M1 model / M2 robust / M3 cluster) ====\n")
ds <- attr(est, "dat.seq")
print(ci_analytic(ds))

cat("\n==== ci_bootstrap (B=15, both) ====\n")
b <- ci_bootstrap(d, B = 15, estimand = "both", seed = 2)
cat("n_fail:", b$n_fail, "\n"); print(b$gamma); print(b$mrd)

cat("\n==== ci_jackknife (hd, small n=50) ====\n")
d50 <- dat_sim(n = 50, seed = 3)
jk <- ci_jackknife(d50, estimand = "hd", max_n = Inf)
print(jk$gamma)

cat("\n==== true_values (small pseudo-RCT n=2e4) ====\n")
tv <- true_values(n = 2e4, seed = 5)
print(round(rbind(surv1 = tv$surv1, surv0 = tv$surv0, mrd = tv$mrd), 4))

cat("\n==== ci_lef (M6 pointwise intervals and simultaneous band) ====\n")
lf <- ci_lef(d, B = 15, seed = 2)
print(lf$mrd)
print(lf$band)

cat("\nSMOKE TEST COMPLETE\n")
