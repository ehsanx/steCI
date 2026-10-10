# THE JENSEN ARTEFACT: how much of the headline shortfall is a METRIC choice, not a variance defect?
#
# THE ISSUE. `se_ratio = mean_se / emp_sd` compares a MEAN OF SEs to an SD OF ESTIMATES. Those are not
# the same functional. Exactly:
#     E[se^2] = mean_se^2 + var(se) = mean_se^2 (1 + CV^2),  CV = sd(se)/mean(se)
#  => rms_se = mean_se * sqrt(1 + CV^2)
#  => se_ratio = mean_se/emp_sd = (rms_se/emp_sd) * 1/sqrt(1 + CV^2)
# The second factor is PURE JENSEN: it drags se_ratio below the truth whenever the SE varies across
# replicates, with NO variance defect required at all.
#
# WHICH RATIO IS RIGHT? An UNBIASED variance estimator satisfies E[se_hat^2] = Var(est), i.e.
# **rms_se/emp_sd = 1**. That -- not mean_se/emp_sd -- is the variance-calibration measure. mean_se is a
# convenient summary, not the estimand.
#
# CAVEAT (stated, not hidden): coverage depends on the JOINT law of (est, se_hat), not on either margin.
# So rms_se/emp_sd is the right VARIANCE-UNBIASEDNESS measure, but coverage remains the ground truth and
# is reported alongside. A previously measured Spearman corr(se_i, |est_i - mean_est|) of
# 0.016/0.063/0.101/0.119/0.105 by tau is weakly POSITIVE, i.e. mildly HELPFUL to coverage.
#
# Guards match simulation/metrics.R (mrd: |est| <= 1, se <= 0.5), so this is comparable to the promoted
# metrics_ci.csv.
#
#   Rscript dev/jensen_decomposition.R
suppressWarnings(suppressMessages({library(dplyr); library(tidyr)}))
ROOT <- "."
# RESULTS PATH. This script reads data/results. The three aggregations:
#   data/results          <- THIS: 81 scenarios, fixed code (shipped)
#   (a discarded 3-scenario, n_sim=4 smoke stub)  <- NOT shipped
#   simulation/output/results                     <- OLD pre-fix run, SUPERSEDED (regenerate locally)
# This script previously read the OLD path. Pooled quantities agree to <=0.001, but cell-level values
# move by up to 0.13 in se_ratio, so the canonical path is mandatory.
RES <- Sys.getenv("RESULTS_DIR",
                  file.path(ROOT, "data/results"))
if (!file.exists(file.path(RES, "metrics_ci.csv"))) stop("no metrics_ci.csv in ", RES)
cat("[jensen] results dir:", RES, "\n")
source(file.path(ROOT, "simulation/scenarios.R"))
rp  <- readRDS(file.path(RES, "replicates_all.rds"))
ci  <- read.csv(file.path(RES, "metrics_ci.csv"), stringsAsFactors = FALSE)
scn <- make_scenarios()[, c("scenario", "n", "gamma.L", "gamma.0", "alpha.0")]

cat("######## THE JENSEN DECOMPOSITION ########\n\n")

# emp_sd from the point rows, guarded exactly as metrics.R does.
pt <- rp %>% filter(method == "point", estimand == "mrd", is.finite(estimate), abs(estimate) <= 1) %>%
  group_by(scenario, tau) %>% summarise(emp_sd = sd(estimate), n_pt = n(), .groups = "drop")

# SE moments per (scenario, tau, method), same guard.
se <- rp %>% filter(method == "m6_analytic_cluster", estimand == "mrd",
                    is.finite(se), se > 0, se <= 0.5) %>%
  group_by(scenario, tau) %>%
  summarise(mean_se = mean(se), rms_se = sqrt(mean(se^2)), cv_se = sd(se) / mean(se),
            n_se = n(), .groups = "drop")

d <- pt %>% inner_join(se, by = c("scenario", "tau")) %>% inner_join(scn, by = "scenario") %>%
  mutate(se_ratio_published = mean_se / emp_sd,          # what the paper's metric computes
         se_ratio_rms       = rms_se / emp_sd,           # Jensen-free: E[se^2] vs Var(est)
         jensen_factor      = 1 / sqrt(1 + cv_se^2))     # the multiplicative artefact

## ---- V. validation: reproduce the published se_ratio ----
v <- ci %>% filter(estimand == "mrd", method == "m6_analytic_cluster") %>%
  select(scenario, tau, se_ratio_pub = se_ratio) %>% inner_join(d, by = c("scenario", "tau"))
cat("V. VALIDATION: my mean_se/emp_sd vs the promoted metrics_ci.csv se_ratio\n")
cat("   rows:", nrow(v), "| max abs diff =", format(max(abs(v$se_ratio_pub - v$se_ratio_published), na.rm = TRUE), digits = 3), "\n")
cat("   identity check  se_ratio_published == se_ratio_rms * jensen_factor : max err =",
    format(max(abs(d$se_ratio_published - d$se_ratio_rms * d$jensen_factor), na.rm = TRUE), digits = 3), "\n\n")

## ---- 1. the headline cell ----
cat("1. THE HEADLINE (gamma.L=3, tau=5, mean over the 9 alpha.0 x gamma.0 cells), by n\n")
h <- d %>% filter(gamma.L == 3, tau == 5) %>% group_by(n) %>%
  summarise(cv_se = mean(cv_se), jensen_factor = mean(jensen_factor),
            published = mean(se_ratio_published), jensen_free = mean(se_ratio_rms), .groups = "drop") %>%
  mutate(jensen_share_of_log_shortfall = -log(jensen_factor) / -log(published))
print(as.data.frame(h %>% mutate(across(where(is.numeric), ~round(., 4)))))
cat("\n   Read: 'published' is what the paper reports. 'jensen_free' is the variance-calibration measure.\n")
cat("   'jensen_share' is the fraction of the log-shortfall that is PURELY the mean-vs-rms metric choice.\n")

## ---- 2. does the n-TREND survive? (our strongest surviving novel claim) ----
cat("\n2. *** DOES THE n-TREND SURVIVE THE CORRECTION? *** (the claim Zhou is a foil for)\n")
tr <- h %>% select(n, published, jensen_free)
cat("   published  :", paste(sprintf("%.3f", tr$published), collapse = " -> "),
    sprintf("  (total drop %.3f)", tr$published[1] - tr$published[nrow(tr)]), "\n")
cat("   jensen-free:", paste(sprintf("%.3f", tr$jensen_free), collapse = " -> "),
    sprintf("  (total drop %.3f)", tr$jensen_free[1] - tr$jensen_free[nrow(tr)]), "\n")
cat("   => the trend", if (tr$jensen_free[1] > tr$jensen_free[nrow(tr)]) "SURVIVES" else "DIES",
    "; it is", round(100 * (1 - (tr$jensen_free[1] - tr$jensen_free[nrow(tr)]) /
                            (tr$published[1] - tr$published[nrow(tr)])), 1), "% smaller than published\n")
cat("   NB CV(se_hat) GROWS with n, so Jensen eats MORE of the shortfall at large n --\n")
cat("   i.e. part of the apparent WIDENING is the metric, not the estimator.\n")

## ---- 3. by horizon ----
cat("\n3. by horizon (gamma.L=3, n=5000)\n")
bt <- d %>% filter(gamma.L == 3, n == 5000) %>% group_by(tau) %>%
  summarise(cv_se = round(mean(cv_se), 3), jensen = round(mean(jensen_factor), 4),
            published = round(mean(se_ratio_published), 3), jensen_free = round(mean(se_ratio_rms), 3),
            jensen_share = round(mean(-log(jensen_factor)) / mean(-log(se_ratio_published)), 3),
            .groups = "drop")
print(as.data.frame(bt))

## ---- 4. the gate cell ----
cat("\n4. THE GATE CELL (gamma.L=0.5, gamma.0=0 -- benign weights, 13% shortfall) and its twin\n")
g <- d %>% filter(gamma.L == 0.5, tau == 5, n == 5000, alpha.0 == 0.25) %>%
  select(gamma.0, cv_se, jensen_factor, published = se_ratio_published, jensen_free = se_ratio_rms)
print(as.data.frame(g %>% mutate(across(where(is.numeric), ~round(., 4)))))
cat("   Does the gate cell's shortfall survive the Jensen correction? If jensen_free ~ 1.0 there,\n")
cat("   the cell that killed every mechanism was a METRIC artefact all along.\n")

## ---- 5. whole-grid summary ----
cat("\n5. across the whole grid (mrd, m6_analytic_cluster): how much is metric?\n")
w <- d %>% filter(se_ratio_published < 1) %>%
  summarise(cells = n(),
            median_jensen_share = round(median(-log(jensen_factor) / -log(se_ratio_published)), 3),
            q25 = round(quantile(-log(jensen_factor) / -log(se_ratio_published), 0.25), 3),
            q75 = round(quantile(-log(jensen_factor) / -log(se_ratio_published), 0.75), 3),
            cells_fully_explained = sum(se_ratio_rms >= 0.98))
print(as.data.frame(w))

saveRDS(d, file.path(ROOT, "dev/jensen_decomposition_results.rds"))
cat("\nJENSEN DONE (saved dev/jensen_decomposition_results.rds)\n")
