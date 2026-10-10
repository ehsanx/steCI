# RE-AGGREGATE the existing 1000-rep run with the guarded metrics.R. NO ARC, NO re-run.
#
# metrics.R is AGGREGATION, not simulation: it reads the stacked replicate rows, which is
# exactly what replicates_all.rds already holds. So the emp_sd / se_ratio contamination
# found in an earlier analysis is fixable today, without touching the simulation.
#
# Writes to *_guarded.* -- it does NOT clobber metrics_ci.csv / metrics.rds. Compare first,
# then promote deliberately.
#
#   Rscript dev/reaggregate_guarded.R
suppressWarnings(suppressMessages({library(dplyr); library(tidyr)}))
ROOT <- "."
RES  <- file.path(ROOT, "simulation/output/results")
source(file.path(ROOT, "simulation/scenarios.R"))
source(file.path(ROOT, "simulation/metrics.R"))

reps <- readRDS(file.path(RES, "replicates_all.rds"))
tv   <- readRDS(file.path(ROOT, "simulation/output/truth_cache.rds"))
truth_fun <- truth_lookup(tv)
old  <- read.csv(file.path(RES, "metrics_ci.csv"), stringsAsFactors = FALSE)

cat("######## RE-AGGREGATION WITH GUARDED metrics.R ########\n\n")
cat("replicate rows:", nrow(reps), "\n\n")

m <- summarize_metrics(reps, truth_fun)
new <- m$ci

## ---- V. VALIDATION: unguarded columns must reproduce the published file exactly ----
cat("V1. the RAW (unguarded) columns must reproduce metrics_ci.csv bit-for-bit:\n")
j <- old %>% select(scenario, estimand, tau, method, se_ratio_pub = se_ratio,
                    emp_sd_pub = emp_sd, mean_se_pub = mean_se, coverage_pub = coverage) %>%
  inner_join(new %>% select(scenario, estimand, tau, method, se_ratio_raw, emp_sd_raw,
                            mean_se_raw, coverage),
             by = c("scenario", "estimand", "tau", "method"))
cat("   rows joined:", nrow(j), "of", nrow(old), "published\n")
cat("   max |emp_sd_raw  - published| =", format(max(abs(j$emp_sd_raw  - j$emp_sd_pub),  na.rm = TRUE), digits = 3), "\n")
cat("   max |mean_se_raw - published| =", format(max(abs(j$mean_se_raw - j$mean_se_pub), na.rm = TRUE), digits = 3), "\n")
cat("   max |se_ratio_raw- published| =", format(max(abs(j$se_ratio_raw- j$se_ratio_pub),na.rm = TRUE), digits = 3), "\n")
cat("   max |coverage    - published| =", format(max(abs(j$coverage    - j$coverage_pub),na.rm = TRUE), digits = 3),
    "  (coverage is deliberately UNguarded)\n\n")

## ---- V2. joint coverage must reproduce the earlier 0.8807 for the band ----
cat("V2. joint (simultaneous) coverage -- must reproduce the earlier band = 0.8807:\n")
jc <- m$joint %>% filter(estimand == "mrd") %>% group_by(method) %>%
  summarise(joint = round(weighted.mean(joint_coverage, n_curves), 4),
            curves = sum(n_curves), .groups = "drop") %>% arrange(joint)
print(as.data.frame(jc))
cat("\n")

## ---- 1. how much contamination was there? ----
cat("1. invalid replicate outputs found (these were silently entering the scale stats):\n")
# NB: count point-estimate failures from m$empirical, which has ONE row per
# (scenario, estimand, tau). The ci table repeats n_invalid_est on every method row, so
# summing it there multiplies the true count by the number of CI methods (~9).
# Bad SEs ARE per-method, so those are counted from the ci table.
inv_est <- m$empirical %>% group_by(estimand) %>%
  summarise(cells = n(), cells_with_bad_est = sum(n_invalid_est > 0),
            total_bad_est = sum(n_invalid_est), .groups = "drop")
inv_se <- new %>% group_by(estimand) %>%
  summarise(method_cells = n(), cells_with_bad_se = sum(n_invalid_se > 0),
            total_bad_se = sum(n_invalid_se), .groups = "drop")
print(as.data.frame(inv_est))
print(as.data.frame(inv_se))
cat("\n   worst |MRD|>1 offenders (point estimates outside the parameter space):\n")
pt <- reps %>% filter(method == "point", estimand == "mrd", is.finite(estimate), abs(estimate) > 1)
print(as.data.frame(pt %>% count(scenario, n, name = "n_impossible") %>%
                      arrange(desc(n_impossible)) %>% head(6)))
cat("   range of the impossible estimates: [", round(min(pt$estimate), 3), ",",
    round(max(pt$estimate), 3), "]\n\n")

## ---- 2. THE HEADLINE: what changed in se_ratio? ----
cat("2. se_ratio: published (contaminated) vs guarded, mrd, mean over 81 scenarios by tau\n")
cmp <- new %>% filter(estimand == "mrd",
                      method %in% c("m6_analytic_cluster", "m6_multiplier", "boot_percentile",
                                    "boot_normal", "m6_band")) %>%
  group_by(method, tau) %>%
  summarise(published = round(mean(se_ratio_raw, na.rm = TRUE), 3),
            guarded   = round(mean(se_ratio, na.rm = TRUE), 3),
            fold      = round(mean(se_ratio_fold, na.rm = TRUE), 3), .groups = "drop")
print(as.data.frame(cmp %>% pivot_longer(c(published, guarded, fold)) %>%
                      pivot_wider(names_from = tau, values_from = value, names_prefix = "tau")))

cat("\n   the cells the guard actually rescues (largest |published - guarded|):\n")
worst <- new %>% filter(estimand == "mrd", is.finite(se_ratio_raw), is.finite(se_ratio)) %>%
  mutate(delta = abs(se_ratio_raw - se_ratio)) %>% arrange(desc(delta)) %>%
  select(scenario, tau, method, se_ratio_raw, se_ratio, se_ratio_fold,
         n_invalid_est, n_invalid_se) %>% head(8)
print(as.data.frame(worst %>% mutate(across(where(is.numeric), ~signif(., 4)))))

## ---- 3. does the earlier headline survive the guard? ----
cat("\n3. does the guard change the earlier story? (m6_analytic_cluster, mrd)\n")
c3 <- new %>% filter(estimand == "mrd", method == "m6_analytic_cluster") %>%
  left_join(reps %>% distinct(scenario, n, gamma.L), by = "scenario") %>%
  filter(gamma.L == 3.0) %>% group_by(n, tau) %>%
  summarise(published = round(mean(se_ratio_raw, na.rm = TRUE), 3),
            guarded   = round(mean(se_ratio, na.rm = TRUE), 3),
            fold      = round(mean(se_ratio_fold, na.rm = TRUE), 3), .groups = "drop")
print(as.data.frame(c3 %>% filter(tau %in% c(1, 5))))
cat("   (the shortfall WORSENS with n at gamma.L=3. If that survives the guard, it is real.)\n")

saveRDS(m, file.path(RES, "metrics_guarded.rds"))
write.csv(m$ci, file.path(RES, "metrics_ci_guarded.csv"), row.names = FALSE)
write.csv(m$empirical, file.path(RES, "metrics_empirical_guarded.csv"), row.names = FALSE)
write.csv(m$joint, file.path(RES, "metrics_joint_guarded.csv"), row.names = FALSE)
cat("\nWROTE (originals untouched): metrics_guarded.rds, metrics_ci_guarded.csv,\n")
cat("      metrics_empirical_guarded.csv, metrics_joint_guarded.csv\n")
cat("\nREAGGREGATE DONE\n")
