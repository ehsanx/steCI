# Two load-bearing checks that were previously taken ON TRUST.
#
# (1) THE 72% POINT MASS AT w == 1. Asserted as a measured fact, but its PROVENANCE was a
#     conjecture ("trial-baseline rows where the cumulative product is still empty"). It matters
#     because Keogh's truncation takes quantile(w, 0.95) over a distribution with this atom, so
#     at gamma.L=3/gamma.0=-2 the 0.95 cap can land ON the atom (P(w<=1) = 0.9524 straddles 0.95).
#     Settle it: decompose the w==1 rows by rownum and by A.baseline.
#
# (2) MARTINUSSEN & SCHEIKE (2006) CONDITION 5.1(a): sup_t E(Y_i W_i^2 X X X) < infinity -- a
#     FINITE-SECOND-MOMENT-OF-THE-WEIGHT condition. Our whole "the theory's conditions hold; our
#     problem is RATE, not VALIDITY" framing rests on it, and we must NOT pretend the weights have
#     infinite variance (they don't). Verify the ALGEBRA against the actual construction, not the
#     DGM prose: if ipw.stab.cum is a finite product of expit-ratios whose linear predictors are
#     linear in a Gaussian-mixture L, then 1/wt.denom = 1 + exp(-eta) is lognormal-mixture-tailed,
#     every moment is finite, and a product of <=5 of them has every moment finite by Holder.
#
#   Rscript dev/verify_weight_structure.R
suppressWarnings(suppressMessages({library(dplyr)}))
ROOT <- "."
for (f in list.files(file.path(ROOT, "steCI/R"), pattern = "[.]R$", full.names = TRUE)) source(f)

mk <- function(n, gamma.L, gamma.0, alpha.0, seed) {
  d <- dat_sim(n = n, n.visit = 5, gamma.0 = gamma.0, gamma.L = gamma.L, alpha.0 = alpha.0,
               alpha.A = -0.04, alpha.L = 0.015, alpha.U = 0.015, seed = seed)
  ste_prepare(d, n.visit = 5, truncate = 0.95)
}

cat("######## (1) PROVENANCE OF THE POINT MASS AT w == 1 ########\n\n")
ds <- mk(1000, 3.0, -2.0, 0.25, 11)
ds$is1 <- abs(ds$ipw.stab.cum - 1) < 1e-12
cat("rows:", nrow(ds), "| P(w == 1) =", round(mean(ds$is1), 4), "\n")
cat("            (this is the atom Keogh's quantile(w, 0.95) is taken over)\n\n")

cat("Decomposition of the w==1 rows -- is the conjecture right?\n")
tab <- ds %>% group_by(rownum_1 = rownum == 1, treated_arm = A.baseline == 1) %>%
  summarise(rows = n(), n_w_eq_1 = sum(is1), pct_of_group_at_1 = round(100 * mean(is1), 1),
            .groups = "drop")
print(as.data.frame(tab))
cat("\n  share of ALL w==1 rows explained by (rownum==1 OR treated arm):",
    round(100 * sum(ds$is1 & (ds$rownum == 1 | ds$A.baseline == 1)) / sum(ds$is1), 2), "%\n")
cat("  ...and are there w==1 rows those two do NOT explain?",
    sum(ds$is1 & ds$rownum != 1 & ds$A.baseline == 0), "\n")
cat("\n  MECHANISM (from the code, not conjecture):\n")
cat("   - ste_prepare.R:81  ipw.stab.lag = ifelse(rownum == 1, 1, lag(ipw.stab)) -> cumprod starts at 1\n")
cat("   - ste_prepare.R:65,73  wt.denom and wt.num are BOTH forced to 1 when A.baseline == 1\n")
cat("     => ipw.stab == 1 on the whole treated arm => cumprod stays exactly 1.\n")

cat("\n  Does the 0.95 cap land ON the atom in this cell?\n")
cat("    P(w <= 1) =", round(mean(ds$ipw.stab.cum <= 1), 4),
    " | quantile(w, 0.95) =", signif(quantile(ds$ipw.stab.cum, 0.95), 6), "\n")
cat("    (if P(w<=1) straddles 0.95 across replicates the cap flips between 1.000 and >1)\n")

cat("\n\n######## (2) M&S CONDITION 5.1(a): IS E[W^2] FINITE? ########\n\n")
cat("Structure check against the ACTUAL construction:\n")
r <- range(ds$wt.denom[ds$A.baseline == 0 & ds$rownum != 5], na.rm = TRUE)
cat("  wt.denom on control, predicted rows: range [", signif(r[1], 4), ",", signif(r[2], 4), "]",
    " -- an expit, so in (0,1), can approach 0 but never equals it\n")
cat("  min wt.denom over ALL rows:", signif(min(ds$wt.denom, na.rm = TRUE), 4),
    " => max 1/wt.denom =", signif(1 / min(ds$wt.denom, na.rm = TRUE), 4), "\n")
cat("  1/wt.denom == 1 + exp(-eta) exactly? ")
eta <- qlogis(pmin(pmax(ds$wt.denom, 1e-12), 1 - 1e-12))
cat(isTRUE(all.equal(1 / ds$wt.denom, 1 + exp(-eta), tolerance = 1e-8)), "\n")
cat("  distinct factors in the cumulative product (max rownum):", max(ds$rownum, na.rm = TRUE),
    "=> a product of at most 5 terms\n")

cat("\nTHE ALGEBRA (this, not a moment estimate, is what settles 5.1(a)):\n")
cat("  wt.denom = expit(b0 + b1*L.lead1), so 1/wt.denom = 1 + exp(-b0 - b1*L.lead1).\n")
cat("  L is a finite Gaussian MIXTURE (dat_sim.R:46 -- A enters the mean), so exp(-b1*L) is a\n")
cat("  finite mixture of lognormals => EVERY moment finite. wt.num is bounded in (0,1).\n")
cat("  Hence ipw.stab <= 1/wt.denom has every moment finite, and ipw.stab.cum is a product of\n")
cat("  at most 5 such factors => every moment finite by Holder.  => E[W^2] < infinity: 5.1(a) HOLDS.\n")

cat("\nEmpirical corroboration (NOT proof -- an alpha=infinity control mimics these):\n")
for (nn in c(1000, 5000)) {
  d2 <- mk(nn, 3.0, -2.0, 0.25, 7)
  w <- d2$ipw.stab.cum
  cat(sprintf("  n=%5d  E[w]=%.4f  E[w^2]=%.4f  E[w^4]=%.3f  sd(log w)=%.4f  max=%.2f\n",
              nn, mean(w), mean(w^2), mean(w^4), sd(log(w)), max(w)))
}
cat("\n  Read: a stable sd(log w) with a growing max is a SLOW-approach signature, not an\n")
cat("  infinite moment. Condition 5.1(a) is a POPULATION statement and it holds; it says nothing\n")
cat("  about how fast the approximation kicks in. That gap IS our finding.\n")
cat("\nWEIGHT_STRUCTURE DONE\n")
