# Two items, both on the CANONICAL aggregation.
#
# (1) RECOMPUTE "only 4 of 81 scenarios reach nominal simultaneous coverage".
#     The old figure came from the pre-fix run, where 18 of 81 cells were short (worst 640/1000) and
#     replicate loss correlates POSITIVELY with measured coverage (r=+0.26),
#     so three of those four cells were themselves short. Recompute honestly.
#
# (2) PER-METHOD CONSTRUCTION-FAILURE RATES. Promised at the paper and never delivered,
#     because the old blanket tryCatch deleted a replicate for EVERY method whenever M6 failed, so
#     per-method rates were unidentifiable. The fix gives each method its own tryCatch, so the
#     rates are now recoverable: a method's failures show up as a deficit in its own n_sim.
#
#   Rscript dev/failure_rates_and_nominal.R
#
# Item (1) is a coverage, so it depends on the truth: it defaults to the exact-truth aggregation
# data/exact/results, which gives the paper's figures (simultaneous band coverage 0.754-0.956,
# median 0.887, 3 of 81 scenarios at >= 0.95). RESULTS_DIR=data/results reproduces the earlier
# Monte-Carlo-truth figures. Item (2) is truth-free and identical under either.
suppressWarnings(suppressMessages({library(dplyr); library(tidyr)}))
ROOT <- Sys.getenv("REPO_ROOT", ".")
RES  <- Sys.getenv("RESULTS_DIR", file.path(ROOT, "data/exact/results"))
source(file.path(ROOT, "simulation/scenarios.R"))
ci  <- read.csv(file.path(RES, "metrics_ci.csv"), stringsAsFactors = FALSE)
scn <- make_scenarios()[, c("scenario","n","alpha.0","gamma.L","gamma.0","event_rate","prevalence","confounding")]
ci  <- merge(ci, scn, by = "scenario")
cat("SOURCE:", RES, "\n\n")

## ---------------- (1) simultaneous coverage: how many scenarios reach nominal? ----------------
j <- read.csv(file.path(RES, "metrics_joint.csv"), stringsAsFactors = FALSE) %>%
  filter(estimand == "mrd", method == "m6_band") %>% merge(scn, by = "scenario")
cat("=== (1) BAND SIMULTANEOUS COVERAGE ACROSS THE 81 SCENARIOS ===\n")
cat("  scenarios:", nrow(j), "| curves:", sum(j$n_curves), "\n")
cat("  pooled (curve-weighted):", round(weighted.mean(j$joint_coverage, j$n_curves), 4), "\n")
cat("  min / median / max      :", paste(round(quantile(j$joint_coverage, c(0,.5,1)), 3), collapse = " / "), "\n")
cat("  reaching >= 0.95        :", sum(j$joint_coverage >= 0.95), "of", nrow(j), "\n")
cat("  reaching >= 0.94        :", sum(j$joint_coverage >= 0.94), "\n")
cat("  reaching >= 0.90        :", sum(j$joint_coverage >= 0.90), "\n")
cat("\n  Are the high-coverage scenarios the SHORT ones? (the old artefact)\n")
cat("  correlation(joint_coverage, n_curves) =",
    round(cor(j$joint_coverage, j$n_curves), 3), "\n")
top <- j %>% arrange(desc(joint_coverage)) %>% head(6) %>%
  select(scenario, n, joint_coverage, n_curves)
print(as.data.frame(top %>% mutate(joint_coverage = round(joint_coverage, 3))))
cat("\n  by n:\n")
print(as.data.frame(j %>% group_by(n) %>%
  summarise(scenarios = n(), median = round(median(joint_coverage), 3),
            n_at_nominal = sum(joint_coverage >= 0.95), .groups = "drop")))

## ---------------- (2) per-method construction-failure rates ----------------
cat("\n\n=== (2) PER-METHOD CONSTRUCTION-FAILURE RATES ===\n")
cat("A method's failures appear as a deficit in its own n_sim (1000 attempted per cell).\n")
cat("Under the OLD blanket tryCatch every method shared one n_sim, so this was unidentifiable.\n\n")
fr <- ci %>% filter(estimand == "mrd") %>%
  group_by(method) %>%
  summarise(cells = n(), attempted = 1000 * n(), completed = sum(n_sim),
            fail_pct = round(100 * (1 - sum(n_sim) / (1000 * n())), 3),
            worst_cell_n = min(n_sim), .groups = "drop") %>% arrange(desc(fail_pct))
print(as.data.frame(fr))

cat("\n  Same, for the constant hazard difference (E1):\n")
fr1 <- ci %>% filter(estimand == "hd") %>% group_by(method) %>%
  summarise(cells = n(), fail_pct = round(100 * (1 - sum(n_sim) / (1000 * n())), 3),
            worst_cell_n = min(n_sim), .groups = "drop") %>% arrange(desc(fail_pct))
print(as.data.frame(fr1))

cat("\n  WHERE the failures concentrate (mrd, methods with any failure):\n")
bad <- ci %>% filter(estimand == "mrd", n_sim < 1000) %>%
  group_by(scenario, n, event_rate, prevalence, confounding) %>%
  summarise(methods_affected = n_distinct(method),
            min_n_sim = min(n_sim), max_n_sim = max(n_sim), .groups = "drop") %>%
  arrange(min_n_sim)
print(as.data.frame(head(bad, 8)))
cat("\n  scenarios with ANY method below 1000:", nrow(bad), "of 81\n")
cat("  scenarios where methods DIFFER in n_sim (proof the per-method guard is live):",
    sum(bad$min_n_sim != bad$max_n_sim), "\n")
cat("\nDONE\n")
