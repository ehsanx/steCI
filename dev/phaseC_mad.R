# Is the late-horizon SE shortfall TAIL-DRIVEN?
#
# Hypothesis (from phaseC_diagnose.R): at late tau a MINORITY of replicates with an
# unstable accumulated B-hat(tau) (shrinking risk set), pushed through the nonlinear
# exp() map, produce extreme MRD-hat. If so, emp_sd is inflated by that tail and BOTH
# the analytic IF and the bootstrap "under-cover" against a target that is itself a
# tail-driven artefact -- i.e. a late-horizon ESTIMATOR instability, not a
# variance-estimator defect.
#
# Decisive quantity: se_ratio computed against a ROBUST scale instead of emp_sd.
#   mean_se / emp_sd    = 0.797 at tau=5   (known)
#   mean_se / robust_sd = ?                 (if ~1.0 => the SE is right for the BULK)
#
# Zero cost: read-only over replicates_all.rds + metrics_ci.csv. No re-run.
#   Rscript dev/phaseC_mad.R
suppressWarnings(suppressMessages({library(dplyr); library(tidyr)}))
RES <- "./simulation/output/results"
ci  <- read.csv(file.path(RES, "metrics_ci.csv"), stringsAsFactors = FALSE)
rp  <- readRDS(file.path(RES, "replicates_all.rds"))

cat("############ is the tau=5 spread tail-driven? ############\n\n")
cat("methods present:", paste(sort(unique(rp$method)), collapse = ", "), "\n\n")

# metrics.R derives emp_sd from the method=="point" rows ONLY, and every CI method
# is then scored against that same emp_sd. So the point rows are the exact basis.
pt <- rp %>% filter(estimand == "mrd", method == "point", is.finite(estimate))
cat("point rows (mrd):", nrow(pt), "| scenarios:", n_distinct(pt$scenario), "\n\n")

# ---------------------------------------------------------------- scales per cell
safe_div <- function(a, b) ifelse(is.finite(b) & b > 0, a / b, NA_real_)

scales <- pt %>%
  group_by(scenario, n, alpha.0, gamma.0, gamma.L, tau) %>%
  summarise(
    n_rep   = n(),
    emp_sd  = sd(estimate),
    mad_sd  = mad(estimate),               # 1.4826 * median|x - med|, normal-consistent
    iqr_sd  = IQR(estimate) / 1.349,       # normal-consistent
    sd_tr1  = { q <- quantile(estimate, c(0.005, 0.995));
                sd(estimate[estimate >= q[1] & estimate <= q[2]]) },   # 1% trimmed
    sd_tr5  = { q <- quantile(estimate, c(0.025, 0.975));
                sd(estimate[estimate >= q[1] & estimate <= q[2]]) },   # 5% trimmed
    kurt    = { z <- estimate - mean(estimate)
                mean(z^4) / (mean(z^2)^2) - 3 },                       # excess kurtosis
    max_z   = safe_div(max(abs(estimate - median(estimate))), mad(estimate)),
    frac_5mad = { m <- mad(estimate)
                  if (is.finite(m) && m > 0) mean(abs(estimate - median(estimate)) > 5 * m) else NA_real_ },
    # share of the total (about-mean) sum of squares contributed by the worst 1% of reps
    var_share_top1 = { z2 <- (estimate - mean(estimate))^2
                       k <- max(1L, ceiling(0.01 * length(z2)))
                       sum(sort(z2, decreasing = TRUE)[seq_len(k)]) / sum(z2) },
    mean_neg_haz = mean(neg_haz_frac, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(sd_over_mad = safe_div(emp_sd, mad_sd),
         sd_over_iqr = safe_div(emp_sd, iqr_sd))

# ---------------------------------------------------------------- validation first
# Discipline: reproduce a published number before trusting anything.
val <- ci %>% filter(estimand == "mrd", method == "m6_analytic_cluster") %>%
  select(scenario, tau, emp_sd_pub = emp_sd) %>%
  left_join(scales %>% select(scenario, tau, emp_sd), by = c("scenario", "tau"))
cat("V. VALIDATION -- my emp_sd vs metrics_ci.csv emp_sd (must be ~0):\n")
cat("   max abs diff =", format(max(abs(val$emp_sd_pub - val$emp_sd), na.rm = TRUE), digits = 3),
    "| cells compared =", sum(is.finite(val$emp_sd_pub - val$emp_sd)), "\n\n")

# ---------------------------------------------------------------- C2a
cat("C2a. emp_sd vs ROBUST scales, by horizon (mean over 81 scenarios)\n")
cat("     (1.00 = Gaussian bulk; >1 = emp_sd inflated by a heavy tail)\n")
a <- scales %>% group_by(tau) %>%
  summarise(emp_sd = mean(emp_sd), mad_sd = mean(mad_sd), iqr_sd = mean(iqr_sd),
            sd_over_mad = mean(sd_over_mad, na.rm = TRUE),
            sd_over_iqr = mean(sd_over_iqr, na.rm = TRUE),
            excess_kurt = mean(kurt, na.rm = TRUE),
            frac_beyond_5mad = mean(frac_5mad, na.rm = TRUE),
            var_share_top1 = mean(var_share_top1, na.rm = TRUE),
            .groups = "drop")
print(as.data.frame(a %>% mutate(across(where(is.numeric), ~round(., 4)))))

# ---------------------------------------------------------------- C2b  THE TEST
cat("\nC2b. *** THE DECISIVE TEST *** se_ratio against emp_sd vs against robust scales\n")
cat("     m6_analytic_cluster. If mean_se/robust ~ 1.0 while mean_se/emp_sd = 0.80,\n")
cat("     the SE is CORRECT for the bulk and emp_sd is tail-inflated.\n")
mse <- ci %>% filter(estimand == "mrd", method == "m6_analytic_cluster",
                     is.finite(mean_se)) %>%
  select(scenario, tau, mean_se)
b <- scales %>% inner_join(mse, by = c("scenario", "tau")) %>%
  group_by(tau) %>%
  summarise(`mean_se/emp_sd`  = mean(safe_div(mean_se, emp_sd), na.rm = TRUE),
            `mean_se/mad_sd`  = mean(safe_div(mean_se, mad_sd), na.rm = TRUE),
            `mean_se/iqr_sd`  = mean(safe_div(mean_se, iqr_sd), na.rm = TRUE),
            `mean_se/sd_tr1`  = mean(safe_div(mean_se, sd_tr1), na.rm = TRUE),
            `mean_se/sd_tr5`  = mean(safe_div(mean_se, sd_tr5), na.rm = TRUE),
            .groups = "drop")
print(as.data.frame(b %>% mutate(across(where(is.numeric), ~round(., 3)))))

cat("\n     same, split by n (does it hold at every sample size?)\n")
b2 <- scales %>% inner_join(mse, by = c("scenario", "tau")) %>%
  group_by(n, tau) %>%
  summarise(emp = round(mean(safe_div(mean_se, emp_sd), na.rm = TRUE), 3),
            mad = round(mean(safe_div(mean_se, mad_sd), na.rm = TRUE), 3), .groups = "drop") %>%
  pivot_longer(c(emp, mad)) %>%
  pivot_wider(names_from = tau, values_from = value, names_prefix = "tau")
print(as.data.frame(b2))

# ---------------------------------------------------------------- C2c
cat("\nC2c. the bootstrap too (boot_normal): is IT also fine against a robust scale?\n")
mseb <- ci %>% filter(estimand == "mrd", method == "boot_normal", is.finite(mean_se)) %>%
  select(scenario, tau, mean_se_boot = mean_se)
cc <- scales %>% inner_join(mseb, by = c("scenario", "tau")) %>%
  filter(mean_se_boot < 5) %>%     # guard: the degenerate 18.62 cell (see C2e)
  group_by(tau) %>%
  summarise(n_cells = n(),
            `boot/emp_sd` = round(mean(safe_div(mean_se_boot, emp_sd), na.rm = TRUE), 3),
            `boot/mad_sd` = round(mean(safe_div(mean_se_boot, mad_sd), na.rm = TRUE), 3),
            .groups = "drop")
print(as.data.frame(cc))

# ---------------------------------------------------------------- C2d mechanism
cat("\nC2d. MECHANISM: are the extreme replicates the negative-hazard-truncated ones?\n")
pt5 <- pt %>% filter(tau == 5) %>% group_by(scenario) %>%
  mutate(m = median(estimate), s = mad(estimate),
         extreme = is.finite(s) & s > 0 & abs(estimate - m) > 5 * s) %>%
  ungroup()
d <- pt5 %>% group_by(extreme) %>%
  summarise(n_reps = n(),
            mean_neg_haz = round(mean(neg_haz_frac, na.rm = TRUE), 4),
            median_neg_haz = round(median(neg_haz_frac, na.rm = TRUE), 4),
            mean_abs_mrd = round(mean(abs(estimate)), 4),
            max_abs_mrd = round(max(abs(estimate)), 4), .groups = "drop")
print(as.data.frame(d))
ok <- pt5 %>% filter(is.finite(neg_haz_frac))
if (nrow(ok) > 10 && sd(ok$neg_haz_frac) > 0) {
  cat("     corr(|MRD-hat - median|, neg_haz_frac) =",
      round(cor(abs(ok$estimate - ok$m), ok$neg_haz_frac, method = "spearman"), 3),
      "(Spearman)\n")
}

# ---------------------------------------------------------------- C2e data bug
cat("\nC2e. DATA BUG #2 -- which cells emit the degenerate bootstrap SE (mean 18.62)?\n")
bug <- ci %>% filter(estimand == "mrd", method == "boot_normal", tau == 5) %>%
  arrange(desc(mean_se)) %>% select(scenario, tau, mean_se, emp_sd, se_ratio, n_sim) %>% head(5)
print(as.data.frame(bug %>% mutate(across(where(is.numeric), ~signif(., 4)))))
cat("\n     ...and the se_ratio=86.6 family (boot_percentile, worst cells):\n")
bug2 <- ci %>% filter(estimand == "mrd", method == "boot_percentile") %>%
  arrange(desc(se_ratio)) %>% select(scenario, tau, mean_se, emp_sd, se_ratio, n_sim) %>% head(5)
print(as.data.frame(bug2 %>% mutate(across(where(is.numeric), ~signif(., 4)))))
cat("\n     ARE THEY THE SAME CELLS AS THE TAIL-DRIVEN ONES? worst sd/mad cells at tau=5:\n")
worst <- scales %>% filter(tau == 5) %>% arrange(desc(sd_over_mad)) %>%
  select(scenario, n, n_rep, emp_sd, mad_sd, sd_over_mad, kurt, max_z) %>% head(5)
print(as.data.frame(worst %>% mutate(across(where(is.numeric), ~signif(., 4)))))

saveRDS(list(scales = scales, by_tau = a, se_ratio_robust = b, by_n = b2,
             boot = cc, mechanism = d, worst = worst),
        "./dev/phaseC_mad_results.rds")

cat("\nINTERPRETATION KEY:\n")
cat(" * mean_se/mad_sd ~ 1.0 at tau=5  => SE is right for the BULK; emp_sd is tail-inflated\n")
cat("   => reframe as a late-horizon ESTIMATOR instability (scope statement), NOT an IF bug.\n")
cat(" * mean_se/mad_sd still ~0.80     => a genuine variance-estimator shortfall; keep digging\n")
cat("   (and the tail hypothesis is dead).\n")
cat("\nTAIL-SCALE CHECK DONE (saved dev/phaseC_mad_results.rds)\n")
