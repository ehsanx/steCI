# Fraction of truncated (negative) additive-hazard evaluations by baseline hazard alpha_0.
# Backs the paper's event-rate statement: "at alpha_0 = 0.05, 14-43% of hazard evaluations are
# truncated, depending on treatment prevalence and confounding (mean 27%); at most 0.2% at
# alpha_0 >= 0.15". Part 1 simulates the 9 prevalence x confounding cells at alpha_0 = 0.05 and
# 0.15; part 2 (only if the per-replicate record data/results/replicates_all.rds is present --
# 91 MB, not in git, see the README) summarises the main simulation's own per-replicate values.
# The neg_haz_frac attribute is computed inside dat_sim() (steCI/R/dat_sim.R).
#   Rscript dev/truncation_by_alpha0.R      (about a minute; writes nothing but stdout)
# Stored output: dev/truncation_by_alpha0.log.
suppressMessages({library(dplyr); library(survival)})
for (f in list.files("steCI/R", pattern = "[.]R$", full.names = TRUE)) source(f)
n <- 4000L
grid <- expand.grid(gamma.0 = c(0, -1, -2), gamma.L = c(0.5, 1.5, 3.0))
for (a0 in c(0.05, 0.15)) {
  fr <- vapply(seq_len(nrow(grid)), function(i) {
    d <- dat_sim(n = n, gamma.0 = grid$gamma.0[i], gamma.L = grid$gamma.L[i], alpha.0 = a0, seed = 2000 + i)
    attr(d, "neg_haz_frac")
  }, 0)
  cat(sprintf("alpha_0 = %.2f: neg_haz_frac over the 9 prevalence x confounding cells (n = %d each)\n", a0, n))
  cat(sprintf("  mean %.4f  min %.4f  max %.4f\n", mean(fr), min(fr), max(fr)))
  print(data.frame(grid, neg_haz_frac = round(fr, 4)), row.names = FALSE)
  cat("\n")
}

# ---- part 2: the main simulation's per-replicate truncation fractions (alpha_0 = 0.15/0.25/0.40) ----
REP <- Sys.getenv("TRUNC_REPLICATES", "data/results/replicates_all.rds")
if (file.exists(REP)) {
  r <- readRDS(REP)
  p <- r[r$method == "point" & r$estimand == "hd", ]          # one row per (scenario, replicate)
  cat(sprintf("Main simulation (%s, neg_haz_frac per replicate; %d replicates):
", REP, nrow(p)))
  for (a0 in sort(unique(p$alpha.0))) {
    x <- p$neg_haz_frac[p$alpha.0 == a0]
    cat(sprintf("  alpha_0 = %.2f: mean %.2g, max %.3g
", a0, mean(x, na.rm = TRUE), max(x, na.rm = TRUE)))
  }
} else {
  cat(sprintf("(part 2 skipped: %s not found)
", REP))
}
