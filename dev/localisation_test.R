# THE LOCALISATION TEST -- the last surviving structural claim from the five dead mechanisms.
#
# THE CLAIM (from an internal check, NEVER TESTED): "per-component sd ratios at tau=5 are 0.804 (B0) /
# 0.841 (B_A) / 0.986 (B_L), and V_IF/V_delta ~ 1 => the delta step and the standardisation are
# EXONERATED; 100% of the defect is the Aalen sandwich for Cov(B_hat), in the INTERCEPT and TREATMENT
# cumulative coefficients at late t -- NOT in B_L."
#
# WHY IT MATTERS. It is a statement about WHERE the defect is, not WHY -- so it does not rescue HC0
# (refuted) or any weight story (structurally excluded: sum_i h_i = p is weight-invariant). But
# "we localise the defect to Cov(B_hat); the standardisation and the delta step are exonerated" is a
# genuine, publishable NARROWING even with no mechanism. It is the only survivor worth testing.
#
# DESIGN.
#  * B = (B0, B_A, B_L): the Aalen cumulative coefficients (intercept, A.baseline, L.baseline).
#  * Per component j and horizon s: the cluster-robust IF SE is sqrt(sum_i IFc[i, j, s]^2)
#    (the same construction ci_lef.R:98 uses for the MRD).
#  * ratio_j = SE_j / sd(B_j across replicates). ratio ~ 1 => that component's variance is calibrated.
#  * ALSO report the MRD ratio from the SAME fits, so the components and the target are comparable.
#
# TWO LESSONS FROM THIS SESSION ARE BUILT IN:
#  1. JENSEN: mean_se/emp_sd is up to 36% metric artefact. `E[se^2] = mean_se^2 (1+CV^2)`, so the
#     calibration measure is **rms_se/emp_sd** (an unbiased variance estimator has E[se^2]=Var).
#     Both are reported; the rms one is the answer.
#  2. >= 24 SEEDS before any number is reported. Three seeds produced a 23-fold error today.
#
#   Rscript dev/localisation_test.R          (LOC_R overrides reps; default 300)
suppressWarnings(suppressMessages({library(dplyr); library(tidyr); library(parallel)}))
ROOT <- "."
R_REPS <- as.integer(Sys.getenv("LOC_R", "300"))
TAU <- 1:5
CELLS <- data.frame(
  label   = c("strong_conf", "gate_cell"),
  gamma.L = c(3.0, 0.5),
  gamma.0 = c(-2.0, 0.0),
  alpha.0 = c(0.25, 0.25), stringsAsFactors = FALSE)
NS <- c(1000L, 5000L)

# One replicate -> B_hat(tau) per component, its IF-based cluster SE, and the MRD + its SE from the
# SAME fit, so component and target are directly comparable.
one <- function(nn, gl, g0, a0, seed, tau_) {
  d  <- dat_sim(n = nn, n.visit = 5, gamma.0 = g0, gamma.L = gl, alpha.0 = a0,
                alpha.A = -0.04, alpha.L = 0.015, alpha.U = 0.015, seed = seed)
  ds <- ste_prepare(d, n.visit = 5, truncate = 0.95)
  fit <- try(.aalen_if(ds, t.hor = tau_, weight_col = "ipw.stab.cum"), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  L0 <- d$L[d$visit == 1]
  comp <- c("B0", "B_A", "B_L")
  out <- list()
  for (s in seq_along(tau_)) {
    B <- fit$B[s, ]
    # per-component cluster-robust SE: the same rss-of-per-subject-IF construction as ci_lef.R:98
    for (j in 1:3) out[[length(out) + 1L]] <- data.frame(
      quantity = comp[j], tau = tau_[s], est = B[j], se = sqrt(sum(fit$IFc[, j, s]^2)))
    # the MRD from the same fit, via the delta method on the same IFc
    s1 <- mean(exp(-(B[1] + B[2] + B[3] * L0))); s0 <- mean(exp(-(B[1] + B[3] * L0)))
    g <- c(mean(exp(-(B[1] + B[2] + B[3]*L0))) - mean(exp(-(B[1] + B[3]*L0))),
           mean(exp(-(B[1] + B[2] + B[3]*L0))),
           mean(L0 * exp(-(B[1] + B[2] + B[3]*L0))) - mean(L0 * exp(-(B[1] + B[3]*L0))))
    ifm <- fit$IFc[, , s] %*% g
    out[[length(out) + 1L]] <- data.frame(quantity = "MRD", tau = tau_[s],
                                          est = s0 - s1, se = sqrt(sum(ifm^2)))
  }
  res <- do.call(rbind, out); res$seed <- seed; res$n <- nn; res
}

cat("######## LOCALISATION TEST ########\n")
cat("R =", R_REPS, "| cells:", nrow(CELLS), "| n:", paste(NS, collapse = ","), "\n")
nw <- max(1L, parallel::detectCores() - 1L); cat("workers:", nw, "\n\n")
cl <- parallel::makeCluster(nw); on.exit(parallel::stopCluster(cl), add = TRUE)
# Windows PSOCK workers inherit NOTHING: export every global by name, and pass values explicitly
# (a bare `n` in a worker resolves to dplyr::n -- "non-numeric matrix extent").
parallel::clusterExport(cl, c("ROOT", "one", "TAU"), envir = environment())
invisible(parallel::clusterEvalQ(cl, {
  suppressWarnings(suppressMessages({library(dplyr); library(tidyr); library(timereg)
                                     library(survival); library(rlang)}))
  for (f in list.files(file.path(ROOT, "steCI/R"), pattern = "[.]R$", full.names = TRUE)) source(f)
  NULL
}))

all <- list(); k <- 0L
for (ci in seq_len(nrow(CELLS))) for (nn in NS) {
  cc <- CELLS[ci, ]; t0 <- proc.time()[["elapsed"]]
  out <- parallel::parLapply(cl, seq_len(R_REPS), function(s, n_, gl_, g0_, a0_, tau_)
    one(n_, gl_, g0_, a0_, 9e6 + s, tau_),
    n_ = nn, gl_ = cc$gamma.L, g0_ = cc$gamma.0, a0_ = cc$alpha.0, tau_ = TAU)
  out <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
  out$cell <- cc$label
  k <- k + 1L; all[[k]] <- out
  cat(sprintf("  [%s] n=%5d reps=%3d  %.1f s\n", cc$label, nn,
              length(unique(out$seed)), proc.time()[["elapsed"]] - t0)); flush.console()
}
r <- do.call(rbind, all)
saveRDS(r, file.path(ROOT, "dev/localisation_results.rds"))

s <- r %>% group_by(cell, n, quantity, tau) %>%
  summarise(R = n(), emp_sd = sd(est), mean_se = mean(se), rms_se = sqrt(mean(se^2)),
            cv_se = sd(se) / mean(se), .groups = "drop") %>%
  mutate(ratio_published = mean_se / emp_sd,   # the paper's metric (Jensen-contaminated)
         ratio_rms = rms_se / emp_sd)          # THE calibration measure

cat("\n1. *** THE CLAIM *** per-component SE/SD at tau=5 (claimed 0.804 B0 / 0.841 B_A / 0.986 B_L)\n")
cat("   ratio_rms is the Jensen-free measure and is the answer; ratio_published is the paper's metric.\n")
print(as.data.frame(s %>% filter(tau == 5) %>%
  select(cell, n, quantity, R, ratio_published, ratio_rms, cv_se) %>%
  mutate(across(where(is.numeric), ~round(., 3)))))

cat("\n2. by horizon -- does B_L stay clean while B0/B_A degrade? (strong_conf)\n")
print(as.data.frame(s %>% filter(cell == "strong_conf") %>%
  select(n, quantity, tau, ratio_rms) %>%
  pivot_wider(names_from = quantity, values_from = ratio_rms) %>%
  mutate(across(where(is.numeric), ~round(., 3)))))

cat("\n3. THE GATE CELL -- where does ITS defect live? (benign weights, genuine ~12% shortfall)\n")
print(as.data.frame(s %>% filter(cell == "gate_cell") %>%
  select(n, quantity, tau, ratio_rms) %>%
  pivot_wider(names_from = quantity, values_from = ratio_rms) %>%
  mutate(across(where(is.numeric), ~round(., 3)))))

cat("\nINTERPRETATION KEY:\n")
cat(" * B_L ~ 1.0 while B0/B_A << 1  -> localisation CONFIRMED: the defect is in Cov(B_hat), and\n")
cat("   specifically in the intercept/treatment components. The delta step and standardisation are\n")
cat("   exonerated (MRD's ratio should then be explicable from B0/B_A alone).\n")
cat(" * ALL components ~ equally short -> the defect is not localised; the claim dies.\n")
cat(" * MRD ratio MUCH worse than every component -> the defect IS in the delta/standardisation step,\n")
cat("   i.e. the opposite of the claim.\n")
cat("\nLOCALISATION DONE (saved dev/localisation_results.rds)\n")
