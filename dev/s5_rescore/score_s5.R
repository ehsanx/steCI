# Score Web Table S5 (appendix tab:stdvar) against the EXACT truths.
#
# Input : per-replicate estimates persisted by dev/s5_rescore/run_s5_reps.R
#         (dev/s5_rescore/reps/cell<ci>_r<a>-<b>.rds; no truth inside them).
# Truth : dev/truth_cache_exact.rds (closed form).
# Output: dev/s5_rescore/s5_scored.rds, dev/s5_rescore/s5_summary.csv, and a printed report.
#
# Arithmetic path of the printed statistic is reproduced exactly from
# dev/experiment_stdvar.R: cov[method, tau] = sum(hit) / nv over the replicates that
# succeeded, then rowMeans over the 5 horizons, then round(., 3). With 400 replicates
# the horizon-mean coverage is a multiple of 1/2000, so an odd numerator sits exactly on
# a 3-decimal rounding tie; such cells are FLAGGED (the printed digit is then a
# floating-point tie-break, and the 4-decimal value is the honest number).
#
# A REPRODUCTION DIAGNOSTIC (section 5) scores the SAME replicates against the
# SUPERSEDED single-draw MC truth data/truth_cache.rds, only to confirm that
# the replicate stream reproduces the printed table (so every change is attributable to
# the truth swap alone). That file is never used for any reported number.
#
#   Rscript dev/s5_rescore/score_s5.R                # requires the full replicate set
#   S5_PARTIAL=1 Rscript dev/s5_rescore/score_s5.R   # interim look on an incomplete set
suppressMessages(library(stats))
ROOT    <- Sys.getenv("REPO_ROOT", ".")
REPS    <- Sys.getenv("S5_REPS", file.path(ROOT, "dev/s5_rescore/reps"))  # override: tests only
T_EXACT <- file.path(ROOT, "dev/truth_cache_exact.rds")
T_OLD   <- file.path(ROOT, "data/truth_cache.rds")   # diagnostic ONLY (sec. 5)
# The table as first printed (scored against the superseded truth), for the before/after columns;
# the original script parsed it from the paper's appendix source, which is not part of this repository.
PRINTED_BEFORE <- data.frame(n = c(300L, 300L, 300L, 1000L, 1000L, 1000L, 5000L),
                             ev = c("low", "med", "high", "low", "med", "high", "med"),
                             p_fixed  = c(0.922, 0.916, 0.934, 0.943, 0.910, 0.933, 0.939),
                             p_random = c(0.922, 0.916, 0.934, 0.943, 0.910, 0.933, 0.939))
OUTDIR  <- Sys.getenv("S5_OUT", file.path(ROOT, "dev/s5_rescore"))
PARTIAL <- identical(Sys.getenv("S5_PARTIAL"), "1")
TH <- 1:5; z <- stats::qnorm(0.975)
METH <- c("cl_fixed", "cl_random", "mult_fixed", "mult_random")

# design -- IDENTICAL to dev/experiment_stdvar.R lines 23-28 and run_s5_reps.R
cells <- data.frame(
  ci      = 1:7,
  n       = c(300L, 300L, 300L, 1000L, 1000L, 1000L, 5000L),
  alpha.0 = c(0.15, 0.25, 0.40, 0.15,  0.25,  0.40,  0.25),
  reps    = c(400L, 400L, 400L, 400L,  400L,  400L,  200L))
cells$ev <- c("low", "med", "high")[match(cells$alpha.0, c(0.15, 0.25, 0.40))]

# ---- 1. load the persisted replicates and assert completeness ----------------------------
files <- list.files(REPS, pattern = "^cell[0-9]+_r[0-9]{3}-[0-9]{3}[.]rds$", full.names = TRUE)
if (!length(files)) stop("no replicate chunks in ", REPS)
est <- do.call(rbind, lapply(files, readRDS))
stopifnot(all(est$ci %in% cells$ci), all(est$tau %in% TH))
key <- paste(est$ci, est$rep, est$tau)
if (anyDuplicated(key)) stop("duplicate (cell, rep, tau) rows across chunks")
incomplete <- character(0)
for (i in seq_len(nrow(cells))) {
  got <- sort(unique(est$rep[est$ci == cells$ci[i]]))
  if (!identical(as.integer(got), seq_len(cells$reps[i])))
    incomplete <- c(incomplete, sprintf("cell %d: %d/%d reps", i, length(got), cells$reps[i]))
  nt <- table(est$rep[est$ci == cells$ci[i]])
  if (length(nt) && any(nt != length(TH))) stop("cell ", i, ": a replicate lacks some horizons")
  # the chunk's design columns must agree with the design table (catches a mislabelled chunk)
  sub <- est[est$ci == cells$ci[i], ]
  if (nrow(sub)) stopifnot(all(sub$n == cells$n[i]),
                           all(abs(sub$alpha0 - cells$alpha.0[i]) <= 1e-9 * cells$alpha.0[i]),
                           all(sub$seed == 1000L * cells$ci[i] + sub$rep))
}
if (length(incomplete)) {
  msg <- paste("INCOMPLETE replicate set:", paste(incomplete, collapse = "; "))
  if (!PARTIAL) stop(msg, "\n(set S5_PARTIAL=1 for an interim, clearly-labelled look)")
  cat("*** PARTIAL / INTERIM -- NOT A RESULT ***\n", msg, "\n\n", sep = "")
}
cat(sprintf("replicate rows: %d  (status: %s)\n", nrow(est),
            paste(names(table(est$status)), table(est$status), sep = "=", collapse = ", ")))
stopifnot(all(is.finite(as.matrix(est[est$status == "ok",
                                      c("mrd", "se_f", "se_r", "qf_lo", "qf_hi", "qr_lo", "qr_hi")]))))

# ---- 2. truth lookup: numeric key match with magnitude-keyed tolerance (never format()) --
truth_table <- function(path, require_exact) {
  tv <- readRDS(path)
  kn <- suppressWarnings(as.numeric(names(tv)))
  stopifnot(!anyNA(kn))
  out <- lapply(cells$alpha.0, function(a0) {
    k <- which(abs(kn - a0) <= 1e-9 * max(1, abs(a0)))
    if (length(k) != 1L) stop("truth key-space in ", path, " does not cover alpha.0 = ", a0,
                              " uniquely (keys: ", paste(names(tv), collapse = ","), ")")
    x <- tv[[k]]
    if (require_exact && !isTRUE(x$exact)) stop("truth for alpha.0=", a0, " is not flagged exact")
    stopifnot(length(x$t.hor) == length(TH), all(abs(x$t.hor - TH) <= 1e-12),
              length(x$mrd) == length(TH), all(is.finite(x$mrd)))
    x$mrd
  })
  do.call(rbind, out)  # rows = cells, cols = horizons
}
TR_EXACT <- truth_table(T_EXACT, require_exact = TRUE)   # asserted BEFORE any scoring loop
cat("exact truth key-space covers all 7 cells x 5 horizons: OK\n")

# ---- 3. scoring (same inequalities, same arithmetic path as experiment_stdvar.R) ----------
score <- function(TR) {
  ok <- est[est$status == "ok", ]
  ok$truth <- TR[cbind(match(ok$ci, cells$ci), match(ok$tau, TH))]
  tr <- ok$truth
  ok$cl_fixed    <- as.numeric(ok$mrd - z * ok$se_f <= tr & tr <= ok$mrd + z * ok$se_f)
  ok$cl_random   <- as.numeric(ok$mrd - z * ok$se_r <= tr & tr <= ok$mrd + z * ok$se_r)
  ok$mult_fixed  <- as.numeric(ok$qf_lo <= tr & tr <= ok$qf_hi)
  ok$mult_random <- as.numeric(ok$qr_lo <= tr & tr <= ok$qr_hi)
  rows <- lapply(seq_len(nrow(cells)), function(i) {
    s <- ok[ok$ci == cells$ci[i], ]
    nv <- length(unique(s$rep))
    if (nv < 2L) {  # only reachable in S5_PARTIAL mode (completeness is asserted otherwise)
      nm <- c(METH, paste0(METH, "_mcse"), "d_cl", "d_cl_mcse", "d_cl_disc",
              "d_mu", "d_mu_mcse", "d_mu_disc")
      return(data.frame(ci = cells$ci[i], n = cells$n[i], ev = cells$ev[i], nv = nv,
                        t(setNames(rep(NA_real_, length(nm)), nm)), stringsAsFactors = FALSE))
    }
    cov <- matrix(0, length(METH), length(TH), dimnames = list(METH, NULL))
    for (m in METH) for (t in seq_along(TH)) cov[m, t] <- sum(s[[m]][s$tau == TH[t]])
    cov <- cov / nv
    mh <- rowMeans(cov)                                  # printed statistic (unrounded)
    # MC SE of the horizon-mean coverage: replicate-level horizon means (keeps the
    # within-replicate correlation across horizons)
    hb <- sapply(METH, function(m) tapply(s[[m]], s$rep, mean))
    mcse <- apply(hb, 2, sd) / sqrt(nv)
    # paired random - fixed difference (same replicate, same truth)
    dcl <- tapply(s$cl_random - s$cl_fixed, s$rep, mean)
    dmu <- tapply(s$mult_random - s$mult_fixed, s$rep, mean)
    data.frame(ci = cells$ci[i], n = cells$n[i], ev = cells$ev[i], nv = nv,
               t(setNames(mh, METH)),
               t(setNames(mcse, paste0(METH, "_mcse"))),
               d_cl = mean(dcl), d_cl_mcse = sd(dcl) / sqrt(nv),
               d_cl_disc = sum(s$cl_random != s$cl_fixed),
               d_mu = mean(dmu), d_mu_mcse = sd(dmu) / sqrt(nv),
               d_mu_disc = sum(s$mult_random != s$mult_fixed),
               stringsAsFactors = FALSE)
  })
  list(per_rep = ok, cell = do.call(rbind, rows),
       per_tau = do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
         s <- ok[ok$ci == cells$ci[i], ]
         do.call(rbind, lapply(TH, function(t) {
           u <- s[s$tau == t, ]
           data.frame(ci = cells$ci[i], n = cells$n[i], ev = cells$ev[i], tau = t,
                      truth = u$truth[1], mean_mrd = mean(u$mrd), bias = mean(u$mrd) - u$truth[1],
                      cl_fixed = mean(u$cl_fixed), cl_random = mean(u$cl_random),
                      mult_fixed = mean(u$mult_fixed), mult_random = mean(u$mult_random))
         }))
       })))
}
tie_flag <- function(x) { f <- x * 1000 - floor(x * 1000); abs(f - 0.5) < 1e-6 }

S <- score(TR_EXACT)

# ---- 4. printed table (the values as first printed, PRINTED_BEFORE above) --------------------
pr <- cbind(line = NA_integer_, PRINTED_BEFORE)
pr <- pr[match(paste(cells$n, cells$ev), paste(pr$n, pr$ev)), ]
stopifnot(!anyNA(pr$p_fixed))

cmp <- data.frame(line = pr$line, n = cells$n, ev = cells$ev, nv = S$cell$nv,
                  printed_fixed = pr$p_fixed, new_fixed = round(S$cell$cl_fixed, 3),
                  fixed_4dp = sprintf("%.4f", S$cell$cl_fixed),
                  fixed_mcse = round(S$cell$cl_fixed_mcse, 4),
                  printed_random = pr$p_random, new_random = round(S$cell$cl_random, 3),
                  random_4dp = sprintf("%.4f", S$cell$cl_random),
                  tie = ifelse(tie_flag(S$cell$cl_fixed) | tie_flag(S$cell$cl_random), "TIE", ""),
                  stringsAsFactors = FALSE)
cmp$d_fixed  <- round(cmp$new_fixed - cmp$printed_fixed, 3)
cmp$d_random <- round(cmp$new_random - cmp$printed_random, 3)
cat("\n=== WEB TABLE S5 (tab:stdvar) RESCORED AGAINST dev/truth_cache_exact.rds ===\n")
cat("printed_* = the table as first printed; new_* = Web Table S5; coverage = mean over tau = 1..5 of M6 cluster-Wald coverage\n")
print(cmp, row.names = FALSE)
cat("\nPaired random - fixed (horizon-mean; same replicates, same truth):\n")
print(data.frame(n = cells$n, ev = cells$ev,
                 cl_diff = sprintf("%+.4f", S$cell$d_cl), cl_mcse = sprintf("%.4f", S$cell$d_cl_mcse),
                 cl_discordant = S$cell$d_cl_disc,
                 mult_diff = sprintf("%+.4f", S$cell$d_mu), mult_mcse = sprintf("%.4f", S$cell$d_mu_mcse),
                 mult_discordant = S$cell$d_mu_disc,
                 mult_fixed = sprintf("%.4f", S$cell$mult_fixed),
                 mult_random = sprintf("%.4f", S$cell$mult_random)), row.names = FALSE)
cat(sprintf("max |random - fixed| (rounded 3dp columns): %.3f\n",
            max(abs(cmp$new_random - cmp$new_fixed))))
cat(sprintf("max |random - fixed| (unrounded, cluster-Wald): %.4f ; (multiplier): %.4f\n",
            max(abs(S$cell$d_cl)), max(abs(S$cell$d_mu))))
pt <- S$per_tau
cat(sprintf(paste0("per horizon, max over cells of |random - fixed|: cluster-Wald %.4f ; multiplier %.4f",
                   " (at n=%d %s, tau=%d)
"),
            max(abs(pt$cl_random - pt$cl_fixed)), max(abs(pt$mult_random - pt$mult_fixed)),
            pt$n[which.max(abs(pt$mult_random - pt$mult_fixed))], pt$ev[which.max(abs(pt$mult_random - pt$mult_fixed))],
            pt$tau[which.max(abs(pt$mult_random - pt$mult_fixed))]))
cat(sprintf("range of new Fixed column: %.3f - %.3f ; printed: %.3f - %.3f\n",
            min(cmp$new_fixed), max(cmp$new_fixed), min(cmp$printed_fixed), max(cmp$printed_fixed)))

# ---- truth-free diagnostics (must not depend on any truth) ---------------------------------
ok <- est[est$status == "ok", ]
tf <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  s <- ok[ok$ci == cells$ci[i], ]
  per_tau <- lapply(TH, function(t) {
    u <- s[s$tau == t, ]; esd <- sd(u$mrd)
    c(ratio_rf = mean(u$se_r / u$se_f),
      rms_ratio_rf = sqrt(mean(u$se_r^2) / mean(u$se_f^2)),
      cal_f = sqrt(mean(u$se_f^2)) / esd, cal_r = sqrt(mean(u$se_r^2)) / esd,
      w_cl_f = mean(2 * z * u$se_f), w_cl_r = mean(2 * z * u$se_r),
      w_mu_f = mean(u$qf_hi - u$qf_lo), w_mu_r = mean(u$qr_hi - u$qr_lo))
  })
  m <- do.call(rbind, per_tau)
  data.frame(n = cells$n[i], ev = cells$ev[i], t(colMeans(m)),
             max_ratio_rf = max(m[, "ratio_rf"]), stringsAsFactors = FALSE)
}))
cat("\nTruth-free diagnostics (means over tau; rms_se/emp_sd calibration):\n")
print(format(tf, digits = 4), row.names = FALSE)

# ---- 5. REPRODUCTION DIAGNOSTIC against the SUPERSEDED MC truth (never a result) ---------
repro <- NULL
if (file.exists(T_OLD)) {
  TR_OLD <- truth_table(T_OLD, require_exact = FALSE)
  So <- score(TR_OLD)
  repro <- data.frame(n = cells$n, ev = cells$ev,
                      printed_fixed = pr$p_fixed, superseded_fixed = round(So$cell$cl_fixed, 3),
                      printed_random = pr$p_random, superseded_random = round(So$cell$cl_random, 3),
                      fixed_4dp = sprintf("%.4f", So$cell$cl_fixed),
                      random_4dp = sprintf("%.4f", So$cell$cl_random))
  repro$match <- repro$printed_fixed == repro$superseded_fixed &
                 repro$printed_random == repro$superseded_random
  cat("\n=== REPRODUCTION DIAGNOSTIC (superseded single-draw MC truth; NOT a result) ===\n")
  print(repro, row.names = FALSE)
  cat(sprintf("printed table reproduced from the persisted replicates + superseded truth: %d/%d cells\n",
              sum(repro$match), nrow(repro)))
  cat("truth shift exact - superseded (MRD), by cell x tau:\n")
  print(round(TR_EXACT - TR_OLD, 5))
}

# ---- 5b. independent cross-check: canonical main-simulation aggregation (exact truths) ----
# Same DGM cells (gamma.0=-1 'med' prevalence, gamma.L=1.5 'moderate' confounding), different
# seeds/replicates: M6 cluster-Wald (m6_analytic_cluster) and multiplier MRD coverage, mean
# over tau. The rescored Fixed column should sit within joint MC error of these.
CANON <- file.path(ROOT, "data/exact/results/metrics_ci.csv")
xchk <- NULL
if (file.exists(CANON)) {
  mc <- read.csv(CANON, stringsAsFactors = FALSE)
  mc <- mc[mc$estimand == "mrd" & mc$method %in% c("m6_analytic_cluster", "m6_multiplier"), ]
  xchk <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    sc <- sprintf("n%d_ev%s_cfmoderate_pvmed", cells$n[i], cells$ev[i])
    g <- function(meth) {
      u <- mc[mc$scenario == sc & mc$method == meth, ]
      stopifnot(nrow(u) == length(TH), all(abs(sort(u$tau) - TH) <= 1e-12))
      u <- u[order(u$tau), ]
      # the canonical truth must be the exact truth (tolerance keyed to magnitude)
      stopifnot(all(abs(u$truth - TR_EXACT[i, ]) <= 1e-8 * pmax(abs(TR_EXACT[i, ]), 1e-3)))
      c(cov = mean(u$coverage), nsim = u$n_sim[1])
    }
    a <- g("m6_analytic_cluster"); b <- g("m6_multiplier")
    data.frame(n = cells$n[i], ev = cells$ev[i], canon_nsim = a[["nsim"]],
               canon_cl = sprintf("%.4f", a[["cov"]]), s5_fixed = sprintf("%.4f", S$cell$cl_fixed[i]),
               s5_fixed_mcse = sprintf("%.4f", S$cell$cl_fixed_mcse[i]),
               canon_mult = sprintf("%.4f", b[["cov"]]), s5_mult_fixed = sprintf("%.4f", S$cell$mult_fixed[i]),
               stringsAsFactors = FALSE)
  }))
  cat("\n=== CROSS-CHECK vs canonical main-simulation aggregation (exact truths; other seeds) ===\n")
  cat("(canonical truth column asserted equal to dev/truth_cache_exact.rds)\n")
  print(xchk, row.names = FALSE)
}

# ---- 6. persist ----------------------------------------------------------------------------
if (!PARTIAL) {
  saveRDS(list(cells = cells, compare = cmp, cell = S$cell, per_tau = S$per_tau,
               per_rep = S$per_rep, truth_free = tf, reproduction = repro, canon_xcheck = xchk,
               truth_file = T_EXACT, created = format(Sys.time()),
               R = R.version.string), file.path(OUTDIR, "s5_scored.rds"))
  write.csv(cmp, file.path(OUTDIR, "s5_summary.csv"), row.names = FALSE)
  cat("\nwrote dev/s5_rescore/s5_scored.rds and dev/s5_rescore/s5_summary.csv\n")
}
cat("SCORE_S5 DONE\n")
