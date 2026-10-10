# STEP 5: extract every number the manuscript rewrite needs, and verify the code fixes
# actually worked in production. Read-only.
#
#   Rscript dev/phaseD_step5_numbers.R
suppressWarnings(suppressMessages({library(dplyr); library(tidyr)}))
ROOT <- "."
NEW  <- file.path(ROOT, "data/results")
OLD  <- file.path(ROOT, "simulation/output/results")

ci  <- read.csv(file.path(NEW, "metrics_ci.csv"), stringsAsFactors = FALSE)
emp <- read.csv(file.path(NEW, "metrics_empirical.csv"), stringsAsFactors = FALSE)
jt  <- read.csv(file.path(NEW, "metrics_joint.csv"), stringsAsFactors = FALSE)
source(file.path(ROOT, "simulation/scenarios.R"))
scn <- make_scenarios()[, c("scenario", "n", "event_rate", "prevalence", "confounding", "gamma.L", "alpha.0", "gamma.0")]
ci  <- merge(ci, scn, by = "scenario")

cat("################ STEP 5 NUMBERS ################\n\n")
cat("scenarios:", length(unique(ci$scenario)), "| rows:", nrow(ci), "\n")
cat("columns present:", paste(intersect(c("se_ratio","se_ratio_fold","se_ratio_raw","fold_sd","median_se",
    "n_invalid_est","n_invalid_se"), names(ci)), collapse = ", "), "\n\n")

## ---- A. DID THE CODE FIXES WORK IN PRODUCTION? ----
cat("A1. REPLICATE COUNTS -- the deletion bug (run_sim.R blanket tryCatch) is fixed if cells are FULL.\n")
cat("    OLD run: 18 of 81 cells were short (worst 640/1000).\n")
r <- ci %>% filter(estimand == "mrd", tau == 5) %>% group_by(scenario, n) %>%
  summarise(n_sim = max(n_sim), .groups = "drop")
cat("    NEW run: cells =", nrow(r), "| min n_sim =", min(r$n_sim), "| cells < 1000:", sum(r$n_sim < 1000), "\n")
if (any(r$n_sim < 1000)) print(as.data.frame(r %>% filter(n_sim < 1000) %>% arrange(n_sim) %>% head(12)))

cat("\nA2. PER-METHOD survival -- the fix means one method's failure no longer deletes the others.\n")
pm <- ci %>% filter(estimand == "mrd", tau == 5) %>% group_by(scenario) %>%
  summarise(n_methods = n_distinct(method), spread = max(n_sim) - min(n_sim), .groups = "drop")
cat("    scenarios where methods have DIFFERING n_sim (proof the per-method guard is live):",
    sum(pm$spread > 0), "of", nrow(pm), "\n")
cat("    (OLD run: 0 of 81 -- all methods died together. >0 here = the fix is working.)\n")

cat("\nA3. GUARDED-OUT invalid outputs (the |MRD|<=1 / se<=0.5 guards):\n")
if ("n_invalid_est" %in% names(ci))
  cat("    invalid estimates:", sum(unique(ci[, c("scenario","tau","n_invalid_est")])$n_invalid_est),
      "| invalid SEs:", sum(ci$n_invalid_se), "\n")

## ---- B. THE BAND: simultaneous coverage (replaces the 0.8807 from the OLD run) ----
cat("\nB. *** THE BAND *** joint (simultaneous) coverage, weighted by curves\n")
b <- jt %>% filter(estimand == "mrd") %>% group_by(method) %>%
  summarise(joint = round(weighted.mean(joint_coverage, n_curves), 4),
            curves = sum(n_curves), .groups = "drop") %>% arrange(desc(joint))
print(as.data.frame(b))
cat("    OLD run band = 0.8807. If the new band is similar, the finding replicates under fixed code.\n")

## ---- C. E1 table (tab:e1) ----
cat("\nC. E1 (hd): coverage / se_ratio by method x n -- tab:e1 replacement\n")
e1 <- ci %>% filter(estimand == "hd", method %in% c("model","robust","cluster","boot_percentile")) %>%
  group_by(method, n) %>%
  summarise(cov = round(mean(coverage), 3), ser = round(mean(se_ratio), 2),
            ser_fold = round(mean(se_ratio_fold, na.rm = TRUE), 2), .groups = "drop")
print(as.data.frame(pivot_wider(e1, names_from = n, values_from = c(cov, ser, ser_fold))))

## ---- D. E2 table (tab:e2) ----
cat("\nD. E2 (mrd): coverage by method x n -- tab:e2 replacement\n")
e2 <- ci %>% filter(estimand == "mrd") %>% group_by(method, n) %>%
  summarise(cov = round(mean(coverage), 3), .groups = "drop")
print(as.data.frame(pivot_wider(e2, names_from = n, values_from = cov)))

## ---- E. THE JENSEN SPLIT on the new run ----
cat("\nE. JENSEN: published mean_se/emp_sd vs the Jensen-free measure (m6_analytic_cluster, mrd)\n")
if (all(c("se_ratio","se_ratio_fold") %in% names(ci))) {
  j <- ci %>% filter(estimand == "mrd", method == "m6_analytic_cluster") %>%
    group_by(n, tau) %>%
    summarise(published = round(mean(se_ratio, na.rm = TRUE), 3),
              fold = round(mean(se_ratio_fold, na.rm = TRUE), 3), .groups = "drop")
  print(as.data.frame(pivot_wider(j, names_from = tau, values_from = c(published, fold))))
}

## ---- F. horizon decline ----
cat("\nF. coverage by horizon (the 'declines mildly' claim), mrd\n")
h <- ci %>% filter(estimand == "mrd",
                   method %in% c("m6_analytic_cluster","m6_multiplier","m6_band","boot_percentile")) %>%
  group_by(method, tau) %>% summarise(cov = round(mean(coverage), 3), .groups = "drop")
print(as.data.frame(pivot_wider(h, names_from = tau, values_from = cov)))

## ---- G. the gamma.L contrast (the n-trend claim, Jensen-corrected) ----
cat("\nG. se_ratio_fold at tau=5 by gamma.L x n -- the n-trend (Zhou contrast) on the NEW run\n")
g <- ci %>% filter(estimand == "mrd", method == "m6_analytic_cluster", tau == 5) %>%
  group_by(gamma.L, n) %>% summarise(fold = round(mean(se_ratio_fold, na.rm = TRUE), 3), .groups = "drop")
print(as.data.frame(pivot_wider(g, names_from = n, values_from = fold)))
cat("    OLD run at gamma.L=3: 0.801 / 0.722 / 0.669 (WIDENS with n).\n")

cat("\nSTEP5_NUMBERS DONE\n")
