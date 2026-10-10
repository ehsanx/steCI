# PHASE C step 1: locate the ~12% M6 variance shortfall. Zero cost -- metrics_ci.csv
# already stores mean_se, emp_sd and se_ratio PER HORIZON. Read-only.
#   Rscript dev/phaseC_diagnose.R
suppressWarnings(suppressMessages({library(dplyr); library(tidyr)}))
ci <- read.csv("./simulation/output/results/metrics_ci.csv",
               stringsAsFactors = FALSE)
rp <- readRDS("./simulation/output/results/replicates_all.rds")
scn <- rp %>% distinct(scenario, n); rm(rp); invisible(gc())

M <- c("m6_analytic_cluster","m6_multiplier","m6_analytic_row","boot_normal","boot_percentile")
d <- ci %>% filter(estimand == "mrd", method %in% M) %>% left_join(scn, by = "scenario")

cat("############ PHASE C: is the SE shortfall horizon-dependent? ############\n\n")

cat("C1. SE ratio (mean estimated SE / empirical SD) BY HORIZON, mean over 81 scenarios\n")
cat("    (1.00 = calibrated; <1 = SE too small)\n")
r1 <- d %>% filter(is.finite(se_ratio), se_ratio < 5) %>%   # guard: degenerate n=300 cell
  group_by(method, tau) %>% summarise(se_ratio = round(mean(se_ratio), 3), .groups = "drop")
print(as.data.frame(pivot_wider(r1, names_from = tau, values_from = se_ratio, names_prefix = "tau")))

cat("\nC2. same, split by n (m6_analytic_cluster only)\n")
r2 <- d %>% filter(method == "m6_analytic_cluster", is.finite(se_ratio), se_ratio < 5) %>%
  group_by(n, tau) %>% summarise(se_ratio = round(mean(se_ratio), 3), .groups = "drop")
print(as.data.frame(pivot_wider(r2, names_from = tau, values_from = se_ratio, names_prefix = "tau")))

cat("\nC3. mean SE by horizon: M6 analytic vs bootstrap vs the empirical SD\n")
r3 <- d %>% filter(method %in% c("m6_analytic_cluster","boot_normal")) %>%
  group_by(method, tau) %>% summarise(mean_se = mean(mean_se), .groups = "drop") %>%
  pivot_wider(names_from = method, values_from = mean_se)
sd_by_tau <- d %>% filter(method == "m6_analytic_cluster") %>%
  group_by(tau) %>% summarise(emp_sd = mean(emp_sd), .groups = "drop")
r3 <- r3 %>% left_join(sd_by_tau, by = "tau") %>%
  mutate(M6_over_SD   = round(m6_analytic_cluster / emp_sd, 3),
         boot_over_SD = round(boot_normal / emp_sd, 3),
         M6_over_boot = round(m6_analytic_cluster / boot_normal, 3))
print(as.data.frame(r3 %>% mutate(across(where(is.numeric), ~round(., 4)))))

cat("\nINTERPRETATION KEY:\n")
cat(" * If M6_over_SD falls with tau but boot_over_SD stays flat -> the analytic IF\n")
cat("   progressively understates variance at late horizons (D(t)^-1 / shrinking risk set).\n")
cat(" * If BOTH fall -> the point estimator's spread is not what either SE targets.\n")
cat(" * Band note: sup-t is studentised by the SAME SE it scales, so a UNIFORM SE error\n")
cat("   largely cancels in c_alpha; the band inherits the SE shortfall directly.\n")
cat("\nPHASE_C1 DONE\n")
