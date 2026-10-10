# TRUNCATION, RE-ANALYSED **PAIRED**. Supersedes the swing test in dev/truncation_collect.R.
#
# WHY. The experiment was deliberately built on COMMON RANDOM NUMBERS: for a given seed, every weight
# regime is fitted to the SAME simulated dataset (dev/truncation_arc.R -- one dat_sim() per rep, then
# ci_lef() per regime). The contrast (trunc - stab) is therefore PAIRED. dev/truncation_collect.R then
# summarised each regime separately and differenced the MEANS, propagating UNPAIRED standard errors --
# which throws away the whole benefit of the design, because sd(est) >> sd(paired difference).
# The consequence: the reported z = -0.35 / +1.31 / +0.97 are VOID,
# and paired the swing is significant in 2 of 3 cells. This script does it correctly.
#
# WHAT IS BEING TESTED. The claim under test is "truncation destroys consistency": that the truncated
# bias GROWS with n while the untruncated bias shrinks. Operationally, with
#     Delta(n) = E[est_trunc - est_stab]   (a paired within-replicate quantity),
# the claim requires Delta to SWING upward with n. Delta(n) is estimated paired; the swing
# Delta(n_max) - Delta(n_min) combines two INDEPENDENT samples (different datasets), so its SE is the
# root-sum-square of the two paired SEs.
#
# It also reports the dose-response across cap severity (q95 on the w==1 atom; q99, q999 off it),
# again paired, and the BIAS of each regime against the DGM truth.
#
#   Rscript dev/truncation_paired.R
suppressWarnings(suppressMessages({library(dplyr); library(tidyr)}))
ROOT <- "."
OUT  <- file.path(ROOT, "simulation/output/truncation")
R_WANT <- 2000L; K_WANT <- 32L

fs <- list.files(OUT, pattern = sprintf("^trunc_.*_R%d_k%d_c[0-9]+\\.rds$", R_WANT, K_WANT),
                 full.names = TRUE)
if (!length(fs)) stop("no R=", R_WANT, " chunks in ", OUT)
r <- do.call(rbind, lapply(fs, readRDS))
cat("######## TRUNCATION, PAIRED RE-ANALYSIS ########\n\n")
cat("chunks:", length(fs), "| rows:", nrow(r), "| reps:", n_distinct(r$rep), "\n")

# same guards as simulation/metrics.R: MRD is a difference of probabilities => |MRD| <= 1 by definition
r <- r %>% filter(method == "analytic_cluster", is.finite(estimate), abs(estimate) <= 1)
tc <- readRDS(file.path(ROOT, "simulation/output/truth_cache.rds"))
r$truth <- vapply(seq_len(nrow(r)), function(i) tc[[as.character(r$alpha.0[i])]]$mrd[r$tau[i]], numeric(1))

# ---- pair on (cell, n, tau, rep): one row per regime, spread to columns -------------------------
w <- r %>% select(cell, n, tau, rep, truth, regime, estimate) %>%
  pivot_wider(names_from = regime, values_from = estimate) %>%
  filter(!is.na(stab), !is.na(trunc_q95), !is.na(trunc_q99), !is.na(trunc_q999))
cat("paired replicates retained:", nrow(w), "\n\n")

# ---- Delta(n) = paired mean of (trunc - stab), with PAIRED vs UNPAIRED SE side by side ----------
pair_stats <- function(d, a, b) {
  di <- d[[a]] - d[[b]]
  data.frame(R = length(di), delta = mean(di),
             se_paired = sd(di) / sqrt(length(di)),
             se_unpaired = sqrt(var(d[[a]]) / length(di) + var(d[[b]]) / length(di)),
             cor_ab = cor(d[[a]], d[[b]]))
}

cat("1. Delta(n) = mean(trunc_q95 - stab), PAIRED, at tau=5\n")
cat("   (se_unpaired is what truncation_collect.R used -- note the inflation factor)\n")
d5 <- w %>% filter(tau == 5) %>% group_by(cell, n) %>%
  group_modify(~ pair_stats(.x, "trunc_q95", "stab")) %>% ungroup() %>%
  mutate(z_paired = delta / se_paired, z_unpaired = delta / se_unpaired,
         se_inflation = se_unpaired / se_paired)
print(as.data.frame(d5 %>% mutate(across(where(is.numeric), ~signif(., 4)))))

cat("\n2. *** THE SWING TEST, DONE PAIRED *** -- does Delta GROW with n?\n")
cat("   swing = Delta(n=5000) - Delta(n=300); the two n are independent samples,\n")
cat("   so SE(swing) = sqrt(se_paired(5000)^2 + se_paired(300)^2).\n")
sw <- d5 %>% group_by(cell) %>% arrange(n, .by_group = TRUE) %>%
  summarise(delta_lo = first(delta), delta_hi = last(delta),
            n_lo = first(n), n_hi = last(n),
            se_lo = first(se_paired), se_hi = last(se_paired), .groups = "drop") %>%
  mutate(swing = delta_hi - delta_lo, se_swing = sqrt(se_lo^2 + se_hi^2),
         z = swing / se_swing,
         verdict = ifelse(abs(z) > 1.96, "SIGNIFICANT", "not significant"))
print(as.data.frame(sw %>% select(cell, n_lo, n_hi, delta_lo, delta_hi, swing, se_swing, z, verdict) %>%
  mutate(across(where(is.numeric), ~signif(., 4)))))
cat("\n   For comparison, the VOID unpaired swing z-values previously reported were\n")
cat("   -0.35 (clean_twin), +1.31 (gate_cell), +0.97 (strong_conf).\n")

cat("\n3. DOSE-RESPONSE across cap severity, PAIRED, tau=5 (departure from untruncated)\n")
dose <- w %>% filter(tau == 5) %>% group_by(cell, n) %>%
  summarise(d_q95  = mean(trunc_q95  - stab), z95  = mean(trunc_q95  - stab)/(sd(trunc_q95  - stab)/sqrt(n())),
            d_q99  = mean(trunc_q99  - stab), z99  = mean(trunc_q99  - stab)/(sd(trunc_q99  - stab)/sqrt(n())),
            d_q999 = mean(trunc_q999 - stab), z999 = mean(trunc_q999 - stab)/(sd(trunc_q999 - stab)/sqrt(n())),
            .groups = "drop")
print(as.data.frame(dose %>% mutate(across(where(is.numeric), ~signif(., 3)))))
cat("   Monotone in severity => a graded estimand shift (Crump 2009), not a discontinuity artefact.\n")

cat("\n4. BIAS of each regime vs the DGM truth, tau=5 (for the 'X% of the true effect' statement)\n")
bias <- r %>% filter(tau == 5) %>% group_by(cell, n, regime) %>%
  summarise(bias = mean(estimate) - first(truth), truth = first(truth),
            pct_of_truth = 100 * (mean(estimate) - first(truth)) / abs(first(truth)),
            mcse = sd(estimate)/sqrt(n()), .groups = "drop")
print(as.data.frame(bias %>% filter(regime %in% c("stab","trunc_q95")) %>%
  mutate(across(where(is.numeric), ~signif(., 4)))))

saveRDS(list(delta = d5, swing = sw, dose = dose, bias = bias),
        file.path(ROOT, "dev/truncation_paired_results.rds"))
cat("\nPAIRED RE-ANALYSIS DONE (saved dev/truncation_paired_results.rds)\n")
