# PHASE C step 2b: the COVERAGE-RELEVANT scale.
#
# Why this exists. phaseC_mad.R showed emp_sd is tail-inflated (sd/mad = 1.11 at tau=5,
# excess kurtosis 5.7), and that mean_se/mad_sd = 0.872 vs mean_se/emp_sd = 0.797. That
# looks like "the SE is right for the bulk, emp_sd is the artefact".
#
# But MAD is a CENTRAL scale and coverage is a QUANTILE property. With excess kurtosis
# the 2.5/97.5 quantiles sit further out than 1.96 x (central scale), so MAD flatters the
# SE. The coverage-relevant, normal-consistent scale is
#       q_sd = (q97.5 - q2.5) / (2 * 1.96)
# which equals sd for a Gaussian (tau=1 is the calibration check) but tracks the tails
# that actually decide whether an interval covers.
#
# DECISION RULE at tau=5:
#   mean_se/q_sd ~ 0.87 (like MAD)     -> tail reframe SURVIVES; scope statement is honest.
#   mean_se/q_sd ~ 0.80 (like emp_sd)  -> tail reframe is DEAD; the SE is genuinely too
#                                         small where it matters. Keep diagnosing.
#
# Cross-check: a normal-approx CI with ratio r has coverage ~ P(|Z| < 1.96r). Compare that
# prediction under each scale against the OBSERVED coverage. Whichever scale predicts the
# measured coverage is the one that describes reality.
#
#   Rscript dev/phaseC_mad2.R
suppressWarnings(suppressMessages({library(dplyr); library(tidyr)}))
RES <- "./simulation/output/results"
ci  <- read.csv(file.path(RES, "metrics_ci.csv"), stringsAsFactors = FALSE)
rp  <- readRDS(file.path(RES, "replicates_all.rds"))
pt  <- rp %>% filter(estimand == "mrd", method == "point", is.finite(estimate))

cat("######## PHASE C2b: which scale explains the OBSERVED coverage? ########\n\n")

sc <- pt %>% group_by(scenario, n, tau) %>%
  summarise(emp_sd = sd(estimate),
            mad_sd = mad(estimate),
            q_sd   = { q <- quantile(estimate, c(0.025, 0.975)); (q[2] - q[1]) / (2 * qnorm(0.975)) },
            n_rep  = n(), .groups = "drop")

safe <- function(a, b) ifelse(is.finite(b) & b > 0, a / b, NA_real_)

# m6_analytic_cluster: mean_se AND its measured coverage, per cell
m6 <- ci %>% filter(estimand == "mrd", method == "m6_analytic_cluster",
                    is.finite(mean_se), is.finite(coverage)) %>%
  select(scenario, tau, mean_se, coverage, coverage_be)

d <- sc %>% inner_join(m6, by = c("scenario", "tau")) %>%
  mutate(r_emp = safe(mean_se, emp_sd),
         r_mad = safe(mean_se, mad_sd),
         r_q   = safe(mean_se, q_sd))

cat("C2b-1. se_ratio under three scales (mean over 81 scenarios).\n")
cat("       tau=1 is the calibration row: all three MUST be ~1.00 there.\n")
t1 <- d %>% group_by(tau) %>%
  summarise(`mean_se/emp_sd` = round(mean(r_emp, na.rm = TRUE), 3),
            `mean_se/mad_sd` = round(mean(r_mad, na.rm = TRUE), 3),
            `mean_se/q_sd  ` = round(mean(r_q,   na.rm = TRUE), 3), .groups = "drop")
print(as.data.frame(t1))

cat("\nC2b-2. *** THE ADJUDICATION *** predicted coverage P(|Z| < 1.96*r) under each\n")
cat("       scale, vs the OBSERVED coverage. The scale that matches reality wins.\n")
t2 <- d %>% group_by(tau) %>%
  summarise(observed          = round(mean(coverage), 3),
            observed_be       = round(mean(coverage_be), 3),
            pred_from_emp_sd  = round(mean(2 * pnorm(qnorm(0.975) * r_emp) - 1, na.rm = TRUE), 3),
            pred_from_mad_sd  = round(mean(2 * pnorm(qnorm(0.975) * r_mad) - 1, na.rm = TRUE), 3),
            pred_from_q_sd    = round(mean(2 * pnorm(qnorm(0.975) * r_q,  ) - 1, na.rm = TRUE), 3),
            .groups = "drop") %>%
  mutate(err_emp = round(pred_from_emp_sd - observed_be, 3),
         err_mad = round(pred_from_mad_sd - observed_be, 3),
         err_q   = round(pred_from_q_sd   - observed_be, 3))
print(as.data.frame(t2))
cat("\n       (compare against observed_be: bias-eliminated coverage is the like-for-like\n")
cat("        target, since the prediction assumes a correctly centred interval.)\n")

cat("\nC2b-3. how much of the tau=5 shortfall does the TAIL actually explain?\n")
g <- d %>% filter(tau == 5) %>%
  summarise(r_emp = mean(r_emp, na.rm = TRUE), r_mad = mean(r_mad, na.rm = TRUE),
            r_q = mean(r_q, na.rm = TRUE))
cat("       gap against emp_sd  :", round(1 - g$r_emp, 3), "\n")
cat("       gap against q_sd    :", round(1 - g$r_q,   3),
    " <- coverage-relevant residual\n")
cat("       gap against mad_sd  :", round(1 - g$r_mad, 3), " (central scale only)\n")
cat("       share of the emp_sd gap removed by using the quantile scale:",
    round((g$r_q - g$r_emp) / (1 - g$r_emp), 3), "\n")

# Does the M6-vs-bootstrap divergence survive a tail-robust denominator? C1 concluded
# "the bootstrap degrades too => not mainly an IF bug" -- but that used emp_sd. Redo it
# against q_sd, with a HARD guard on the degenerate bootstrap-SE cells (mean_se=1481).
cat("\nC2b-4. M6 vs bootstrap against the quantile scale (C1's conclusion #2 re-tested).\n")
bt <- ci %>% filter(estimand == "mrd", method == "boot_normal", is.finite(mean_se)) %>%
  select(scenario, tau, se_boot = mean_se)
db <- sc %>% inner_join(bt, by = c("scenario", "tau")) %>%
  inner_join(m6 %>% select(scenario, tau, mean_se), by = c("scenario", "tau")) %>%
  # guard: drop cells where the bootstrap SE is absurd relative to the spread
  mutate(bad = !is.finite(se_boot) | se_boot > 5 * emp_sd)
cat("       degenerate bootstrap cells excluded:", sum(db$bad), "of", nrow(db), "\n")
t4 <- db %>% filter(!bad) %>% group_by(tau) %>%
  summarise(n_cells = n(),
            `M6/q_sd`   = round(mean(safe(mean_se, q_sd), na.rm = TRUE), 3),
            `boot/q_sd` = round(mean(safe(se_boot, q_sd), na.rm = TRUE), 3),
            `M6/boot`   = round(mean(safe(mean_se, se_boot), na.rm = TRUE), 3), .groups = "drop")
print(as.data.frame(t4))
cat("\n       If boot/q_sd stays ~1.0 while M6/q_sd decays -> the shortfall IS M6-specific\n")
cat("       (the IF), and C1's 'the bootstrap degrades too' was the tail-inflated denominator.\n")

# The out-of-range estimates: |MRD| > 1 is definitionally impossible.
cat("\nC2b-5. out-of-range point estimates (|MRD-hat| > 1 is impossible by definition)\n")
oor <- pt %>% mutate(bad = abs(estimate) > 1) %>% group_by(tau) %>%
  summarise(n_rep = n(), n_impossible = sum(bad),
            pct = round(100 * mean(bad), 3), .groups = "drop")
print(as.data.frame(oor))
cat("\n       which cells produce them (tau=5)?\n")
oor5 <- pt %>% filter(tau == 5, abs(estimate) > 1) %>%
  count(scenario, n, name = "n_impossible") %>% arrange(desc(n_impossible)) %>% head(8)
print(as.data.frame(oor5))

saveRDS(list(scales = sc, ratios = t1, adjudication = t2, m6_vs_boot = t4, oor = oor,
             oor_cells = oor5), "./dev/phaseC_mad2_results.rds")
cat("\nPHASE_C2b DONE (saved dev/phaseC_mad2_results.rds)\n")
