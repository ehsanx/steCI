# Do q95, q99 and q999 AGREE? -- the check before reporting all three caps together.
#
# The paired re-analysis (dev/truncation_paired.R) ran the swing test only for q95, Keogh's rule.
# Before recommending that the paper report several caps side by side, verify they tell the SAME
# story: if q95 and q99 disagreed on the direction or significance of the n-swing, reporting both
# would create a discrepancy needing explanation rather than a graded trade.
#
# All regimes share the same simulated data (common random numbers), so every contrast is PAIRED.
#
#   Rscript dev/truncation_allcaps.R
suppressWarnings(suppressMessages({library(dplyr); library(tidyr)}))
ROOT <- "."
OUT  <- file.path(ROOT, "simulation/output/truncation")
fs <- list.files(OUT, pattern = "^trunc_.*_R2000_k32_c[0-9]+\\.rds$", full.names = TRUE)
r <- do.call(rbind, lapply(fs, readRDS)) %>%
  filter(method == "analytic_cluster", is.finite(estimate), abs(estimate) <= 1)
tc <- readRDS(file.path(ROOT, "simulation/output/truth_cache.rds"))
r$truth <- vapply(seq_len(nrow(r)), function(i) tc[[as.character(r$alpha.0[i])]]$mrd[r$tau[i]], numeric(1))

w <- r %>% select(cell, n, tau, rep, truth, regime, estimate) %>%
  pivot_wider(names_from = regime, values_from = estimate) %>%
  filter(!is.na(stab), !is.na(trunc_q95), !is.na(trunc_q99), !is.na(trunc_q999))

cat("######## DO THE CAPS AGREE? ########\n\n")
cat("1. PAIRED swing test, run separately for EACH cap (tau=5)\n")
cat("   swing = Delta(n=5000) - Delta(n=300), Delta = paired mean(trunc - stab)\n\n")
out <- list(); k <- 0L
for (cap in c("trunc_q95", "trunc_q99", "trunc_q999")) {
  d <- w %>% filter(tau == 5) %>% group_by(cell, n) %>%
    summarise(delta = mean(.data[[cap]] - stab),
              se = sd(.data[[cap]] - stab) / sqrt(n()), .groups = "drop")
  s <- d %>% group_by(cell) %>% arrange(n, .by_group = TRUE) %>%
    summarise(d_lo = first(delta), d_hi = last(delta),
              se_lo = first(se), se_hi = last(se), .groups = "drop") %>%
    mutate(cap = cap, swing = d_hi - d_lo, se_swing = sqrt(se_lo^2 + se_hi^2),
           z = swing / se_swing, sig = ifelse(abs(z) > 1.96, "SIG", "ns"))
  k <- k + 1L; out[[k]] <- s
}
sw <- do.call(rbind, out)
print(as.data.frame(sw %>% select(cell, cap, d_lo, d_hi, swing, z, sig) %>%
  mutate(across(where(is.numeric), ~signif(., 3)))))

cat("\n2. AGREEMENT CHECK -- same SIGN and same significance verdict across caps, per cell?\n")
ag <- sw %>% group_by(cell) %>%
  summarise(signs = paste(sign(swing), collapse = "/"),
            verdicts = paste(sig, collapse = "/"),
            agree_sign = n_distinct(sign(swing)) == 1,
            agree_verdict = n_distinct(sig) == 1, .groups = "drop")
print(as.data.frame(ag))

cat("\n3. CALIBRATION and COVERAGE by cap (tau=5) -- the numbers the paper would tabulate\n")
cal <- r %>% filter(tau == 5) %>% group_by(cell, n, regime) %>%
  summarise(emp_sd = sd(estimate),
            rms_se = sqrt(mean(se[is.finite(se) & se > 0 & se <= 0.5]^2)),
            coverage = mean(lower <= first(truth) & upper >= first(truth), na.rm = TRUE),
            bias_pct = 100 * (mean(estimate) - first(truth)) / abs(first(truth)),
            .groups = "drop") %>%
  mutate(rms_over_sd = rms_se / emp_sd)
print(as.data.frame(cal %>% filter(n == 5000) %>%
  select(cell, regime, rms_over_sd, coverage, bias_pct) %>%
  mutate(across(where(is.numeric), ~round(., 3)))))

cat("\nREAD: if signs and verdicts agree across caps, the caps tell ONE story of graded severity and\n")
cat("can be tabulated together. If they disagree, reporting several caps needs an explanation.\n")
saveRDS(list(swing = sw, agree = ag, calib = cal), file.path(ROOT, "dev/truncation_allcaps.rds"))
cat("\nALLCAPS DONE\n")
