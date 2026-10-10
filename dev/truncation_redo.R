# TRUNCATION REDO -- production run. Replaces the R=50 evidence behind our sharpest surviving claim.
#
# THE CLAIM: "Truncation restores calibration but destroys consistency: the untruncated bias SHRINKS
# with n (+0.0128/+0.0121/+0.0073 at n=300/1000/5000, tau=4) while the truncated bias GROWS
# (+0.0081/+0.0092/+0.0167 = 17% of the true effect at n=5000)."
#
# WHY IT MUST BE REDONE (dev/probe_truncation_power.R, measured):
#  * R=50 gives the difference-in-differences **power 0.067** -- z = 0.46. A coin flip. The claim is
#    almost certainly noise as it stands.
#  * The cap rule is KEOGH'S (analysis_seq.R:182, not ours): quantile(w, 0.95) taken over a
#    distribution with a ~65% POINT MASS AT EXACTLY 1. Measured: P(w==1) = 0.654, P(w<=1) = 0.959,
#    so the cap is EXACTLY 1.000 -- it deletes ALL up-weighting while keeping ALL down-weighting, and
#    flips across replicates. The "grows with n" result may be an artefact of that discontinuity
#    rather than a property of truncation.
#  * Pairing (common random numbers) buys only 21% (paired SD 0.153 vs 0.193), so R must do the work:
#    **R = 2000 => power 0.83**.
#
# DESIGN. Same data per replicate, five weight regimes, the REAL estimator path (ste_fit ->
# timereg::aalen + mrd_from_aalen). Regimes:
#   stab   = ipw.stab.cum            (untruncated; what the simulation actually uses)
#   keogh  = ipw.stab.cum.trunc      (quantile(w, 0.95) -- the reference rule, lands on the atom)
#   q99    = pmin(w, quantile(w, .99))          -- a plain upper quantile, off the atom
#   tail95 = pmin(w, quantile(w[w>1], .95))     -- UPPER-TAIL-ONLY: the proposed fix
#   abs10  = pmin(w, 10)                        -- a fixed absolute cap
# If "grows with n" holds for keogh but NOT for tail95/abs10, the result is the DISCONTINUITY, not
# truncation. If it holds for all, it is truncation (and Crump 2009 predicts it).
#
#   Rscript dev/truncation_redo.R          (writes dev/truncation_redo_results.rds)
suppressWarnings(suppressMessages({library(dplyr); library(parallel)}))
ROOT <- "."
# REPS ARE ALLOCATED BY COST, NOT FLAT. Measured wall cost on 11 workers: 0.16 / 0.52 / 8.08 s/rep at
# n = 300/1000/5000 -- a 50x spread. n=300 is both the CHEAPEST and the NOISIEST endpoint (paired SD
# falls with n), so it should carry more of the precision. Flat R=2000 costs ~10 h wall and wastes it
# at n=5000; this allocation roughly EQUALISES MCSE(Delta) across the endpoints at ~40% of the cost.
# Target: SE(swing) ~ 0.0045 against a claimed swing of 0.0141 => z ~ 3.1, power ~0.88.
R_BY_N <- c("300" = 3000L, "1000" = 2000L, "5000" = 1000L)
R_SCALE <- as.numeric(Sys.getenv("TRUNC_SCALE", "1"))   # set <1 for a smoke test
TAU <- 1:5
CELLS <- data.frame(gamma.L = c(3.0, 0.5), gamma.0 = c(-2.0, -2.0), alpha.0 = c(0.25, 0.25),
                    label = c("strong_conf", "weak_conf_control"), stringsAsFactors = FALSE)
NS <- c(300L, 1000L, 5000L)

one <- function(n, gamma.L, gamma.0, alpha.0, seed, TAU) {
  d  <- dat_sim(n = n, n.visit = 5, gamma.0 = gamma.0, gamma.L = gamma.L, alpha.0 = alpha.0,
                alpha.A = -0.04, alpha.L = 0.015, alpha.U = 0.015, seed = seed)
  ds <- ste_prepare(d, n.visit = 5, truncate = 0.95)
  w  <- ds$ipw.stab.cum
  ds$w_q99 <- pmin(w, stats::quantile(w, 0.99))
  up <- w[w > 1]
  ds$w_tail95 <- if (length(up) > 20) pmin(w, stats::quantile(up, 0.95)) else w
  ds$w_abs10 <- pmin(w, 10)
  cols <- c(stab = "ipw.stab.cum", keogh = "ipw.stab.cum.trunc",
            q99 = "w_q99", tail95 = "w_tail95", abs10 = "w_abs10")
  res <- lapply(names(cols), function(k) {
    fit <- try(ste_fit(ds, d, t.hor = TAU, weight_col = cols[[k]], estimand = "mrd"), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    data.frame(regime = k, tau = TAU, est = fit$mrd)
  })
  res <- do.call(rbind, res[!vapply(res, is.null, logical(1))])
  if (is.null(res)) return(NULL)
  res$seed <- seed; res$n <- n; res$gamma.L <- gamma.L; res$gamma.0 <- gamma.0
  res$alpha.0 <- alpha.0
  # cap diagnostics: is the cap sitting on the atom this replicate?
  res$cap_keogh <- unname(stats::quantile(w, 0.95))
  res$p_eq_1 <- mean(abs(w - 1) < 1e-12)
  res
}

cat("######## TRUNCATION REDO -- production ########\n")
cat("reps by n:", paste(sprintf("n=%s:R=%d", names(R_BY_N), round(R_BY_N * R_SCALE)), collapse = "  "),
    "| cells:", nrow(CELLS), "\n")
nw <- max(1L, parallel::detectCores() - 1L)
cat("workers:", nw, "\n\n")
cl <- parallel::makeCluster(nw)
on.exit(parallel::stopCluster(cl), add = TRUE)
# Windows PSOCK workers inherit NOTHING from the master: every global the worker touches -- including
# the function `one` itself -- must be exported by name, or you get
# "11 nodes produced errors; first error: could not find function \"one\"".
parallel::clusterExport(cl, c("ROOT", "one", "TAU"), envir = environment())
invisible(parallel::clusterEvalQ(cl, {
  suppressWarnings(suppressMessages({library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang)}))
  for (f in list.files(file.path(ROOT, "steCI/R"), pattern = "[.]R$", full.names = TRUE)) source(f)
  NULL
}))

all <- list(); k <- 0L
for (ci in seq_len(nrow(CELLS))) for (nn in NS) {
  cc <- CELLS[ci, ]
  t0 <- proc.time()[["elapsed"]]
  R_REPS <- max(2L, as.integer(round(R_BY_N[[as.character(nn)]] * R_SCALE)))
  seeds <- (ci * 1e6L) + (match(nn, NS) * 1e5L) + seq_len(R_REPS)
  # Pass EVERY value explicitly. Do NOT rely on lexical scoping into the worker: `n` there resolves to
  # dplyr::n (a function), which yields "non-numeric matrix extent" from dat_sim's matrix(nrow = n, ...).
  out <- parallel::parLapply(cl, seeds, function(s, n_, gl_, g0_, a0_, tau_)
    one(n_, gl_, g0_, a0_, s, tau_),
    n_ = nn, gl_ = cc$gamma.L, g0_ = cc$gamma.0, a0_ = cc$alpha.0, tau_ = TAU)
  out <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
  el <- proc.time()[["elapsed"]] - t0
  k <- k + 1L; all[[k]] <- out
  cat(sprintf("  [%s] n=%5d  reps done=%4d  elapsed=%6.1f s (%.2f s/rep wall)\n",
              cc$label, nn, length(unique(out$seed)), el, el / R_REPS))
  flush.console()
}
res <- do.call(rbind, all)
saveRDS(res, file.path(ROOT, "dev/truncation_redo_results.rds"))
cat("\nSaved dev/truncation_redo_results.rds :", nrow(res), "rows\n")
cat("TRUNCATION_REDO DONE\n")
