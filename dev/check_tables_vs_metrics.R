# Check the HAND-TYPED manuscript tables against the promoted metrics files.
#
# WHY. tab:e1 and tab:e2 are typed by hand (the auto-generated \input fragments idea was never
# built), so they can drift silently from the results. This script recomputes each printed cell
# from the authoritative exact-truth aggregation data/exact/results/metrics_ci.csv and diffs it
# against the values printed in the paper (Tables 3 and 5 of the revised paper).
#
#   Rscript dev/check_tables_vs_metrics.R      # expect "MISMATCHES ... 0" for both tables
#
# Env: RESULTS_DIR (default data/exact/results).
suppressWarnings(suppressMessages({library(dplyr); library(tidyr)}))
ROOT <- Sys.getenv("REPO_ROOT", ".")
RES  <- Sys.getenv("RESULTS_DIR", file.path(ROOT, "data/exact/results"))
source(file.path(ROOT, "simulation/scenarios.R"))
ci  <- read.csv(file.path(RES, "metrics_ci.csv"), stringsAsFactors = FALSE)
scn <- make_scenarios()[, c("scenario", "n")]
ci  <- merge(ci, scn, by = "scenario")

cat("######## HAND-TYPED TABLES vs PROMOTED metrics_ci.csv ########\n\n")

## ---------------- tab:e1 (the paper) -- E1, coverage + SE ratio ----------------
printed_e1 <- tribble(
  ~method,           ~n,    ~cov,   ~ser,
  "model",            300L, 0.822,  0.70,
  "model",           1000L, 0.785,  0.65,
  "model",           5000L, 0.711,  0.61,
  "robust",           300L, 0.819,  0.68,
  "robust",          1000L, 0.790,  0.63,
  "robust",          5000L, 0.724,  0.60,
  "cluster",          300L, 0.909,  0.89,
  "cluster",         1000L, 0.895,  0.83,
  "cluster",         5000L, 0.855,  0.81,
  "boot_percentile",  300L, 0.921,  0.92,
  "boot_percentile", 1000L, 0.904,  0.85,
  "boot_percentile", 5000L, 0.857,  0.82)

actual_e1 <- ci %>% filter(estimand == "hd", method %in% unique(printed_e1$method)) %>%
  group_by(method, n) %>%
  summarise(cov_actual = mean(coverage), ser_actual = mean(se_ratio), .groups = "drop")

e1 <- printed_e1 %>% left_join(actual_e1, by = c("method", "n")) %>%
  mutate(cov_diff = round(cov - round(cov_actual, 3), 4),
         ser_diff = round(ser - round(ser_actual, 2), 4))
cat("tab:e1 -- E1 constant hazard difference (printed vs recomputed, rounded as printed)\n")
print(as.data.frame(e1 %>% mutate(cov_actual = round(cov_actual, 3), ser_actual = round(ser_actual, 3))))
cat("\n  MISMATCHES (coverage):", sum(e1$cov_diff != 0), " (SE ratio):", sum(e1$ser_diff != 0), "\n\n")

## ---------------- tab:e2 (the paper) -- E2, coverage averaged over tau ----------
printed_e2 <- tribble(
  ~method,                  ~n,    ~cov,
  "m6_analytic_row",         300L, 0.868, "m6_analytic_row",        1000L, 0.859, "m6_analytic_row",        5000L, 0.862,
  "m6_analytic_cluster",     300L, 0.920, "m6_analytic_cluster",    1000L, 0.917, "m6_analytic_cluster",    5000L, 0.921,
  "m6_multiplier",           300L, 0.906, "m6_multiplier",          1000L, 0.902, "m6_multiplier",          5000L, 0.902,
  "m6_band",                 300L, 0.961, "m6_band",                1000L, 0.959, "m6_band",                5000L, 0.960,
  "m6a2_analytic_cluster",   300L, 0.914, "m6a2_analytic_cluster",  1000L, 0.912, "m6a2_analytic_cluster",  5000L, 0.918,
  "m6a2_multiplier",         300L, 0.902, "m6a2_multiplier",        1000L, 0.898, "m6a2_multiplier",        5000L, 0.901,
  "boot_percentile",         300L, 0.929, "boot_percentile",        1000L, 0.925, "boot_percentile",        5000L, 0.923)

actual_e2 <- ci %>% filter(estimand == "mrd", method %in% unique(printed_e2$method)) %>%
  group_by(method, n) %>% summarise(cov_actual = mean(coverage), .groups = "drop")

e2 <- printed_e2 %>% left_join(actual_e2, by = c("method", "n")) %>%
  mutate(cov_diff = round(cov - round(cov_actual, 3), 4))
cat("tab:e2 -- E2 marginal risk difference, coverage averaged over horizons\n")
print(as.data.frame(e2 %>% mutate(cov_actual = round(cov_actual, 3))))
cat("\n  MISMATCHES:", sum(e2$cov_diff != 0), "\n\n")

## ---------------- tab:e1_eventrate (the paper's Table 4) -- E1 coverage by event rate x n ----------
printed_e1ev <- tribble(
  ~event_rate, ~n,    ~model, ~robust, ~cluster, ~boot_percentile,
  "low",   300L, 0.781, 0.781, 0.887, 0.904,
  "low",  1000L, 0.734, 0.742, 0.865, 0.877,
  "low",  5000L, 0.604, 0.617, 0.775, 0.784,
  "med",   300L, 0.823, 0.821, 0.910, 0.922,
  "med",  1000L, 0.793, 0.797, 0.904, 0.914,
  "med",  5000L, 0.723, 0.736, 0.872, 0.874,
  "high",  300L, 0.862, 0.856, 0.929, 0.936,
  "high", 1000L, 0.829, 0.831, 0.916, 0.923,
  "high", 5000L, 0.806, 0.818, 0.918, 0.914) %>%
  pivot_longer(c(model, robust, cluster, boot_percentile), names_to = "method", values_to = "cov")
ev  <- make_scenarios()[, c("scenario", "event_rate")]
actual_e1ev <- merge(ci, ev, by = "scenario") %>%
  filter(estimand == "hd", method %in% unique(printed_e1ev$method)) %>%
  group_by(event_rate, n, method) %>% summarise(cov_actual = mean(coverage), .groups = "drop")
e1ev <- printed_e1ev %>% left_join(actual_e1ev, by = c("event_rate", "n", "method")) %>%
  mutate(cov_diff = round(cov - round(cov_actual, 3), 4))
cat("tab:e1_eventrate -- E1 coverage by event rate and n\n")
print(as.data.frame(e1ev %>% mutate(cov_actual = round(cov_actual, 3))))
cat("\n  MISMATCHES:", sum(e1ev$cov_diff != 0), "\n\n")

## ---------------- what tab:e2 does NOT report, and now must ----------------
cat("NOT IN ANY TABLE -- the E2 SE ratio (this is why the 143.6 contamination never surfaced):\n")
gap <- ci %>% filter(estimand == "mrd",
                     method %in% c("m6_analytic_cluster","m6_multiplier","boot_percentile")) %>%
  group_by(method, n) %>%
  summarise(se_ratio_guarded = round(mean(se_ratio, na.rm = TRUE), 3),
            se_ratio_fold    = round(mean(se_ratio_fold, na.rm = TRUE), 3),
            se_ratio_raw_pub = round(mean(se_ratio_raw, na.rm = TRUE), 2), .groups = "drop")
print(as.data.frame(gap))
cat("\nJoint (simultaneous) coverage -- also absent from every table:\n")
j <- read.csv(file.path(RES, "metrics_joint.csv"), stringsAsFactors = FALSE)
print(as.data.frame(j %>% filter(estimand == "mrd") %>% group_by(method) %>%
  summarise(joint = round(weighted.mean(joint_coverage, n_curves), 3), .groups = "drop") %>%
  arrange(desc(joint))))
## ---------------- Web Tables S3 and S8 (E2 by horizon), recomputed ----------------
s38 <- ci %>% filter(estimand == "mrd",
                     method %in% c("m6_analytic_row", "m6_analytic_cluster", "m6_multiplier", "m6_band",
                                   "m6a2_analytic_cluster", "m6a2_multiplier",
                                   "boot_percentile", "boot_normal", "boot_basic"))
cat("\nWeb Table S3 -- E2 coverage by method and horizon (mean over the 81 scenarios):\n")
print(as.data.frame(s38 %>% group_by(method, tau) %>%
  summarise(v = round(mean(coverage), 3), .groups = "drop") %>%
  pivot_wider(names_from = tau, values_from = v, names_prefix = "tau")))
cat("\nWeb Table S8 -- E2 SE ratio by method and horizon (mean over the 81 scenarios):\n")
print(as.data.frame(s38 %>% group_by(method, tau) %>%
  summarise(v = round(mean(se_ratio, na.rm = TRUE), 2), .groups = "drop") %>%
  pivot_wider(names_from = tau, values_from = v, names_prefix = "tau")))
cat("\nTABLE_CHECK DONE\n")
