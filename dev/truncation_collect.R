# Collect + ANALYSE the truncation experiment. The analysis below is PRE-REGISTERED: it is written
# before the data exist so it cannot be fitted post hoc. Run after the ARC array finishes.
#
#   Rscript dev/truncation_collect.R
suppressWarnings(suppressMessages({library(dplyr); library(tidyr)}))
geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
getc <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) x else d }
ROOT <- getc("REPO_ROOT", ".")
OUT  <- getc("TRUNC_OUT", file.path(ROOT, "simulation/output/truncation"))

# Filenames carry R and CHUNK, so a smoke run cannot be silently mixed into a production collection.
# Default to the production settings; override with TRUNC_R / TRUNC_CHUNK to collect something else.
R_WANT <- geti("TRUNC_R", 2000L); K_WANT <- geti("TRUNC_CHUNK", 32L)
pat <- sprintf("^trunc_.*_R%d_k%d_c[0-9]+\\.rds$", R_WANT, K_WANT)
fs <- list.files(OUT, pattern = pat, full.names = TRUE)
allf <- list.files(OUT, pattern = "^trunc_.*\\.rds$")
if (!length(fs)) stop("no output matching R=", R_WANT, ", CHUNK=", K_WANT, " in ", OUT,
                      "\n  (", length(allf), " other trunc_*.rds present: ",
                      paste(utils::head(allf, 3), collapse = ", "), ")")
if (length(allf) > length(fs))
  cat("NOTE: ignoring", length(allf) - length(fs), "file(s) from other R/CHUNK settings (e.g. smoke runs)\n")
r <- do.call(rbind, lapply(fs, readRDS))
r <- r[!vapply(r, is.null, logical(1))[1] & TRUE, , drop = FALSE]
cat("collected", length(fs), "files |", nrow(r), "rows |", n_distinct(r$rep), "reps\n")
EXPECT <- ceiling(R_WANT / K_WANT) * 9L   # 3 cells x 3 n x n_chunks
if (length(fs) < EXPECT) cat("WARNING: only", length(fs), "of", EXPECT,
                            "expected chunks -- the array may be incomplete. Re-submit to resume.\n")
tc <- readRDS(file.path(ROOT, "simulation/output/truth_cache.rds"))
r$truth <- vapply(seq_len(nrow(r)), function(i) tc[[as.character(r$alpha.0[i])]]$mrd[r$tau[i]], numeric(1))

# Guards match simulation/metrics.R: MRD is a difference of probabilities, so |MRD| <= 1 by
# definition; an SE above 0.5 spans the whole parameter space. Counted, never silently dropped.
r$ok_est <- is.finite(r$estimate) & abs(r$estimate) <= 1
r$ok_se  <- is.finite(r$se) & r$se > 0 & r$se <= 0.5
cat("guarded out:", sum(!r$ok_est), "estimates,", sum(!r$ok_se), "SEs of", nrow(r), "\n\n")

M <- "analytic_cluster"   # the paper's recommended analytic method
d <- r %>% filter(method == M)

fold_sd <- function(e) unname(quantile(abs(e - mean(e)), 0.95) / qnorm(0.975))

# Per (cell, n, regime, tau): bias with MC error, the calibration ratios, coverage.
s <- d %>% filter(ok_est) %>% group_by(cell, gamma.L, gamma.0, n, regime, tau) %>%
  summarise(R = n(), truth = first(truth),
            bias = mean(estimate) - first(truth),
            mcse_bias = sd(estimate) / sqrt(n()),
            emp_sd = sd(estimate), fold = fold_sd(estimate),
            coverage = mean(lower <= first(truth) & upper >= first(truth), na.rm = TRUE),
            mcse_cov = sqrt(mean(lower <= first(truth) & upper >= first(truth), na.rm = TRUE) *
                            (1 - mean(lower <= first(truth) & upper >= first(truth), na.rm = TRUE)) / n()),
            .groups = "drop")
se_s <- d %>% filter(ok_se) %>% group_by(cell, n, regime, tau) %>%
  summarise(mean_se = mean(se), rms_se = sqrt(mean(se^2)), cv_se = sd(se) / mean(se), .groups = "drop")
s <- s %>% left_join(se_s, by = c("cell", "n", "regime", "tau")) %>%
  mutate(se_ratio_published = mean_se / emp_sd,     # the paper's metric
         se_ratio_rms       = rms_se / emp_sd,      # Jensen-free: E[se^2] vs Var(est). THE calibration measure.
         se_ratio_fold      = mean_se / fold)
saveRDS(s, file.path(ROOT, "dev/truncation_summary.rds"))

cat("######## PRE-REGISTERED ANALYSIS ########\n\n")
cat("H1. DOES TRUNCATION RESTORE CALIBRATION?  (se_ratio_rms -> 1.0)\n")
cat("    Jensen-free by construction, so this is a real calibration test, not a metric artefact.\n")
print(as.data.frame(s %>% filter(tau == 5) %>%
  select(cell, n, regime, se_ratio_rms, coverage) %>%
  pivot_wider(names_from = regime, values_from = c(se_ratio_rms, coverage)) %>%
  mutate(across(where(is.numeric), ~round(., 3)))))

cat("\nH2. *** THE CLAIM ***  DOES THE TRUNCATED BIAS GROW WITH n WHILE THE UNTRUNCATED SHRINKS?\n")
cat("    Claimed (R=50, power 0.067): stab 0.0128->0.0073 (shrinks); trunc 0.0081->0.0167 (grows).\n")
b <- s %>% filter(tau %in% c(4, 5)) %>% select(cell, n, regime, tau, bias, mcse_bias, R)
print(as.data.frame(b %>% mutate(across(c(bias, mcse_bias), ~round(., 5)))))

cat("\n    Delta(n) = bias(trunc_q95) - bias(stab), with MC error, by cell and n (tau=5):\n")
dl <- s %>% filter(tau == 5, regime %in% c("stab", "trunc_q95")) %>%
  select(cell, n, regime, bias, mcse_bias) %>%
  pivot_wider(names_from = regime, values_from = c(bias, mcse_bias)) %>%
  mutate(delta = bias_trunc_q95 - bias_stab,
         mcse_delta = sqrt(mcse_bias_trunc_q95^2 + mcse_bias_stab^2),
         z = delta / mcse_delta)
print(as.data.frame(dl %>% mutate(across(where(is.numeric), ~round(., 5)))))
cat("\n    The claim needs Delta to SWING (negative at n=300 -> positive at n=5000). Test the swing:\n")
for (cl in unique(dl$cell)) {
  x <- dl %>% filter(cell == cl) %>% arrange(n)
  if (nrow(x) >= 2) {
    sw <- x$delta[nrow(x)] - x$delta[1]
    se <- sqrt(x$mcse_delta[nrow(x)]^2 + x$mcse_delta[1]^2)
    cat(sprintf("      %-12s swing(n=%d->%d) = %+.5f, SE = %.5f, z = %+.2f  %s\n",
                cl, x$n[1], x$n[nrow(x)], sw, se, sw/se,
                if (abs(sw/se) > 1.96) "SIGNIFICANT" else "NOT SIGNIFICANT"))
  }
}

cat("\nH3. *** IS IT TRUNCATION, OR KEOGH'S CAP LANDING ON THE w==1 ATOM? ***\n")
cat("    trunc_q95 caps at EXACTLY 1.000 (on the atom); q99/q999 cap off it. If the bias appears at\n")
cat("    q95 but NOT at q99/q999, the result is the DISCONTINUITY, not truncation.\n")
print(as.data.frame(s %>% filter(tau == 5) %>% select(cell, n, regime, bias, mcse_bias) %>%
  pivot_wider(names_from = regime, values_from = c(bias, mcse_bias)) %>%
  mutate(across(where(is.numeric), ~round(., 5)))))

cat("\nH4. THE GATE CELL (benign weights, genuine ~12% shortfall): does truncation touch it?\n")
print(as.data.frame(s %>% filter(cell == "gate_cell", tau == 5) %>%
  select(n, regime, bias, se_ratio_rms, coverage) %>%
  mutate(across(where(is.numeric), ~round(., 4)))))

cat("\nH5. the bias-variance trade, made explicit (tau=5)\n")
print(as.data.frame(s %>% filter(tau == 5) %>%
  select(cell, n, regime, bias, emp_sd, coverage) %>%
  mutate(rmse = sqrt(bias^2 + emp_sd^2)) %>%
  mutate(across(where(is.numeric), ~round(., 4)))))

cat("\nCOLLECT DONE (saved dev/truncation_summary.rds)\n")
