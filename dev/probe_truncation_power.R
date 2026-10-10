# PILOT + POWER for the truncation redo. MEASURE FIRST, DESIGN SECOND.
#
# WHY THIS EXISTS. The rate-run plan died because it was designed without a power calculation and
# then found to have power 0.074. Same mistake is not being made twice: this script measures the
# noise, then computes the R the design actually needs, BEFORE any production run.
#
# THE CLAIM UNDER TEST (currently rests on R=50 -- our SHARPEST surviving novel claim):
#   "Truncation restores calibration but destroys consistency: the UNTRUNCATED bias SHRINKS with n
#    (+0.0128/+0.0121/+0.0073 at n=300/1000/5000, tau=4) while the TRUNCATED bias GROWS
#    (+0.0081/+0.0092/+0.0167 = 17% of the true effect at n=5000)."
#
# TWO KNOWN PROBLEMS WITH THE EXISTING EVIDENCE:
#  1. R = 50. Three seeds produced a 23-fold error elsewhere in this project today.
#  2. The cap rule is Keogh's `quantile(w, 0.95)` (analysis_seq.R:182 -- NOT our invention), taken over
#     a distribution with a ~67% POINT MASS AT EXACTLY 1. Measured at gamma.L=3/gamma.0=-2:
#     P(w<=1) = 0.9507 and quantile(w, 0.95) = EXACTLY 1.000. So the cap deletes ALL up-weighting while
#     keeping ALL down-weighting -- and it FLIPS between 1.000 and >1 across replicates. The
#     "destroys consistency" result may be an artefact of that discontinuity rather than a property of
#     truncation.
#
# THE DESIGN INSIGHT: use COMMON RANDOM NUMBERS. All weight regimes share dat_sim + ste_prepare, so the
# per-replicate difference (est_trunc - est_stab) is PAIRED and far less noisy than either bias alone.
# The regime contrast at each n is therefore cheap; only the ACROSS-n comparison is not.
#
#   Rscript dev/probe_truncation_power.R
suppressWarnings(suppressMessages({library(dplyr)}))
ROOT <- "."
for (f in list.files(file.path(ROOT, "steCI/R"), pattern = "[.]R$", full.names = TRUE)) source(f)
tc <- readRDS(file.path(ROOT, "simulation/output/truth_cache.rds"))

cat("######## TRUNCATION REDO: PILOT + POWER (measure before designing) ########\n\n")
cat("cores available:", parallel::detectCores(), "\n\n")

TAU <- 1:5
# One replicate: shared data, then the SAME estimator (ste_fit -> timereg::aalen + mrd_from_aalen,
# the REAL pipeline path) under several weight columns. Common random numbers => paired contrasts.
one <- function(n, gamma.L, gamma.0, alpha.0, seed) {
  d  <- dat_sim(n = n, n.visit = 5, gamma.0 = gamma.0, gamma.L = gamma.L, alpha.0 = alpha.0,
                alpha.A = -0.04, alpha.L = 0.015, alpha.U = 0.015, seed = seed)
  ds <- ste_prepare(d, n.visit = 5, truncate = 0.95)   # makes ipw.stab.cum + ipw.stab.cum.trunc
  w  <- ds$ipw.stab.cum
  # alternative cap rules, all applied to the SAME data
  ds$w_q99    <- pmin(w, stats::quantile(w, 0.99))                       # a plain upper quantile
  up <- w[w > 1]
  ds$w_tail95 <- if (length(up) > 20) pmin(w, stats::quantile(up, 0.95)) else w  # UPPER-TAIL-ONLY rule
  ds$w_abs10  <- pmin(w, 10)                                             # fixed absolute cap
  cols <- c(stab = "ipw.stab.cum", keogh = "ipw.stab.cum.trunc",
            q99 = "w_q99", tail95 = "w_tail95", abs10 = "w_abs10")
  out <- lapply(names(cols), function(k) {
    fit <- try(ste_fit(ds, d, t.hor = TAU, weight_col = cols[[k]], estimand = "mrd"), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    data.frame(regime = k, tau = TAU, est = fit$mrd, seed = seed, n = n)
  })
  out <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
  attr(out, "cap") <- c(keogh = unname(stats::quantile(w, 0.95)),
                        q99 = unname(stats::quantile(w, 0.99)),
                        tail95 = if (length(up) > 20) unname(stats::quantile(up, 0.95)) else NA,
                        p_le_1 = mean(w <= 1), p_eq_1 = mean(abs(w - 1) < 1e-12))
  out
}

## ---- pilot ----
NSEED <- 30
cat("PILOT:", NSEED, "seeds at n=1000, gamma.L=3, gamma.0=-2, alpha.0=0.25\n")
t0 <- proc.time()[["elapsed"]]
pil <- do.call(rbind, lapply(1:NSEED, function(s) one(1000L, 3.0, -2.0, 0.25, 1000 + s)))
el <- proc.time()[["elapsed"]] - t0
cat("  elapsed", round(el, 1), "s =>", round(el / NSEED, 2), "s per rep for ALL 5 regimes\n\n")

cap1 <- attr(one(1000L, 3.0, -2.0, 0.25, 1001), "cap")
cat("CAP DIAGNOSTIC (one replicate):\n")
cat("  P(w == 1) =", round(cap1[["p_eq_1"]], 4), " | P(w <= 1) =", round(cap1[["p_le_1"]], 4), "\n")
cat("  Keogh cap = quantile(w, 0.95) =", signif(cap1[["keogh"]], 6),
    if (abs(cap1[["keogh"]] - 1) < 1e-9) "  <-- LANDS EXACTLY ON THE ATOM" else "", "\n")
cat("  q99 cap =", signif(cap1[["q99"]], 6), " | upper-tail-95 cap =", signif(cap1[["tail95"]], 6), "\n\n")

truth5 <- tc[["0.25"]]$mrd[5]
cat("BIAS AND PAIRED CONTRASTS at tau=5 (truth =", round(truth5, 5), "), 30 seeds:\n")
w5 <- pil %>% filter(tau == 5) %>% select(seed, regime, est) %>% tidyr::pivot_wider(names_from = regime, values_from = est)
s <- data.frame(
  regime = c("stab", "keogh", "q99", "tail95", "abs10"),
  bias = sapply(c("stab","keogh","q99","tail95","abs10"), function(k) mean(w5[[k]]) - truth5),
  sd_est = sapply(c("stab","keogh","q99","tail95","abs10"), function(k) sd(w5[[k]])))
s$mcse_bias_R500 <- round(s$sd_est / sqrt(500), 5)
print(as.data.frame(s %>% mutate(across(where(is.numeric), ~signif(., 4)))))

cat("\nPAIRED difference (trunc - stab) -- the payoff of common random numbers:\n")
for (k in c("keogh", "q99", "tail95", "abs10")) {
  dd <- w5[[k]] - w5[["stab"]]
  cat(sprintf("  %-7s mean diff = %+.5f | SD of paired diff = %.5f | MCSE at R=500 = %.5f | ratio vs unpaired = %.2fx\n",
              k, mean(dd), sd(dd), sd(dd)/sqrt(500), sd(dd)/sd(w5[["stab"]])))
}

## ---- power ----
cat("\n######## POWER: what R does the CLAIM actually need? ########\n")
cat("The claim is a DIFFERENCE-IN-DIFFERENCES: does (bias_trunc - bias_stab) GROW with n?\n")
cat("Claimed at tau=4: stab 0.0128 -> 0.0073 (shrinks); keogh 0.0081 -> 0.0167 (grows).\n")
cat("So Delta(n) = bias_trunc - bias_stab goes -0.0047 (n=300) -> +0.0094 (n=5000): a swing of 0.0141.\n\n")
sd_paired <- sd(w5[["keogh"]] - w5[["stab"]])
cat("  measured SD of the paired difference at n=1000, tau=5 =", signif(sd_paired, 4), "\n")
for (R in c(50, 250, 500, 1000, 2000)) {
  se_delta <- sd_paired / sqrt(R)
  se_diff  <- se_delta * sqrt(2)          # two independent n's
  z <- 0.0141 / se_diff
  cat(sprintf("  R=%4d -> MCSE(Delta)=%.5f, SE(swing)=%.5f, z=%.2f, power(5%%)=%.3f%s\n",
              R, se_delta, se_diff, z, round(pnorm(abs(z) - 1.96) , 3),
              if (R == 50) "   <-- WHAT THE CURRENT CLAIM RESTS ON" else ""))
}
cat("\n  NOTE: the paired SD is measured at tau=5/n=1000 and used as a proxy at every n; the production\n")
cat("  run must re-measure per cell. This is a SIZING estimate, not a result.\n")
cat("\nPROBE DONE\n")
