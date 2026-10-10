# Verify the LIVE manuscript claim that the Approach-2 weight-estimation correction changes the
# standard errors "by under 1% throughout the design". This claim is
# FALSE. Check it against the CANONICAL aggregation.
#
#   Rscript dev/check_a2_claim.R
suppressWarnings(suppressMessages({library(dplyr)}))
PHD <- "./data/results"
OLD <- "./simulation/output/results"

run <- function(path, label) {
  ci <- read.csv(file.path(path, "metrics_ci.csv"), stringsAsFactors = FALSE)
  a1 <- ci %>% filter(estimand == "mrd", method == "m6_analytic_cluster") %>%
    select(scenario, tau, se_a1 = mean_se)
  a2 <- ci %>% filter(estimand == "mrd", method == "m6a2_analytic_cluster") %>%
    select(scenario, tau, se_a2 = mean_se)
  d <- inner_join(a1, a2, by = c("scenario", "tau")) %>%
    filter(is.finite(se_a1), is.finite(se_a2), se_a1 > 0) %>%
    mutate(pct = 100 * (se_a2 - se_a1) / se_a1)
  cat("---", label, "(", nrow(d), "cells ) ---\n")
  cat("  median % change :", sprintf("%+.2f", median(d$pct)), "\n")
  cat("  range           : [", sprintf("%+.2f", min(d$pct)), ",", sprintf("%+.2f", max(d$pct)), "]\n")
  cat("  |change| > 1%   :", sprintf("%.1f%%", 100 * mean(abs(d$pct) > 1)),
      sprintf("(%d of %d cells)", sum(abs(d$pct) > 1), nrow(d)), "\n")
  cat("  |change| > 2%   :", sprintf("%.1f%%", 100 * mean(abs(d$pct) > 2)), "\n")
  cat("  by horizon (median %):\n")
  print(as.data.frame(d %>% group_by(tau) %>%
    summarise(median_pct = round(median(pct), 2), max_abs_pct = round(max(abs(pct)), 2),
              frac_over_1pct = round(mean(abs(pct) > 1), 3), .groups = "drop")))
  invisible(d)
}
cat("######## IS 'APPROACH-2 CHANGES SEs BY UNDER 1% THROUGHOUT' TRUE? ########\n\n")
if (dir.exists(PHD)) run(PHD, "CANONICAL results") else cat("(results path not found)\n")
cat("\n")
run(OLD, "OLD pre-fix (superseded; for comparison only)")
cat("\nVERDICT: the claim is true only if |change| never exceeds 1%.\n")
