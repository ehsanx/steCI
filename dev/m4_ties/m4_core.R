# Core of the M4 tie experiment: ONE canonical replicate, bootstrap refitted two ways on the SAME
# resamples of the SAME simulated dataset.
#
#   jitter arm  : the canonical M4 refit, ste_fit() -> timereg::aalen() with read.surv's tie jitter.
#   exact arm   : tie-exact (simultaneous) refit. Primary value = E2 (timereg on the multiplicity-
#                 collapsed resample) where E2 was computed, else E1 (from-scratch simultaneous fit);
#                 E2 is computed on every resample in "full" mode, and in "advance" mode on the first
#                 n_e2_check(n) resamples and on every resample where E1 meets a singular or
#                 ill-conditioned risk-set Gram matrix (so timereg's own singular handling decides).
#
# The replicate is regenerated from the canonical seed stream (run_sim.R: set.seed(20250228);
# sample.int(.Machine$integer.max, 1000); run_one(): set.seed(all_seeds[rep]) -> dat_sim ->
# ste_estimate -> ci_analytic -> ci_bootstrap(seed = NULL, B = 200)), so the resamples are the
# canonical resamples of the main simulation.
#
# MODES
#   full    : runs the canonical call sequence verbatim (ste_estimate, ci_analytic, then the
#             bootstrap loop with the real timereg refits), so the jitter arm is RECOMPUTED and is
#             checked against the canonical replicates_all.rds rows (se, lower, upper) downstream.
#   advance : skips the timereg refits; the RNG is advanced exactly as timereg's read.surv would
#             advance it (runif(nties) once in the hd fit and once in the mrd fit; nothing if the
#             resample's ste_prepare fails, since the fits are never reached). The jitter arm is
#             then read from the canonical replicates_all.rds. Valid only when every canonical
#             refit either succeeded or failed inside ste_prepare -- checked per replicate by
#             requiring (# prepare failures here) == canonical boot_n_fail.
#
# Persisted per replicate (estimates first, truth only at scoring time):
#   draws: B x (gamma, mrd1..5) for jitter (full mode), E1, E2 (NA where not computed);
#   per-resample nties, dt, prepare/fit failure flags, E1 singular counts and min reciprocal
#   condition; the point estimates; mode; timings; the seed. NO truth is read here.

# E2 cross-check budget per replicate in advance mode (E1 = E2 to ~1e-13 wherever E1 is not flagged;
# the check is a running validation, not a source of the primary value). Pilot used 10 at every n.
n_e2_check <- function(n) if (n >= 5000L) 2L else if (n >= 1000L) 3L else 10L
RC_FLAG <- 1e-10

# E1 with timereg's observed singular rule: an exactly singular risk-set Gram matrix contributes a
# ZERO increment (timereg's aalen() returns a flat cumulative coefficient there; see
# validate_tie_exact.R). Flags any singular or ill-conditioned time so the caller can defer to E2.
fit_exact_simul_guarded <- function(ds, dl, t.hor = 1:5, weight_col = "ipw.stab.cum") {
  w <- ds[[weight_col]]
  A <- ds$A.baseline; L <- ds$L.baseline
  start <- ds$time.new; stop <- ds$time.stop.new; ev <- ds$event == 1
  P <- .gram_P(w, A, L)
  et <- sort(unique(stop[ev]))
  fe <- factor(stop[ev], levels = et)
  dN <- rowsum(cbind(w[ev], w[ev] * A[ev], w[ev] * L[ev]), fe, reorder = TRUE)
  # MRD
  G <- .risk_sums(start, stop, P, et)
  sol <- .solve3(G, dN)
  sing <- !is.finite(sol$det) | sol$det == 0
  dB <- sol$z
  r <- dN - cbind(G[, 1] * dB[, 1] + G[, 2] * dB[, 2] + G[, 3] * dB[, 3],
                  G[, 2] * dB[, 1] + G[, 4] * dB[, 2] + G[, 5] * dB[, 3],
                  G[, 3] * dB[, 1] + G[, 5] * dB[, 2] + G[, 6] * dB[, 3])
  dB <- dB + .solve3(G, r)$z
  dB[sing, ] <- 0
  Bc <- apply(dB, 2, cumsum); if (!is.matrix(Bc)) Bc <- matrix(Bc, nrow = 1)
  k <- findInterval(t.hor, et)
  Bh <- matrix(0, length(t.hor), 3); Bh[k > 0, ] <- Bc[k[k > 0], , drop = FALSE]
  L0 <- dl$L[dl$visit == 1]
  surv1 <- vapply(seq_along(t.hor), function(h) mean(exp(-(Bh[h, 1] + Bh[h, 2] + Bh[h, 3] * L0))), 0)
  surv0 <- vapply(seq_along(t.hor), function(h) mean(exp(-(Bh[h, 1] + Bh[h, 3] * L0))), 0)
  rc_m <- sol$rc[!sing]
  # HD (semiparametric), timereg grid {0} U unique(stop)
  g <- sort(unique(c(0, stop[stop > 0]))); gs <- g[-1]; dt <- diff(g)
  G2 <- .risk_sums(start, stop, P, gs)
  m11 <- G2[, 1]; m1L <- G2[, 3]; mLL <- G2[, 6]; m1A <- G2[, 2]; mAL <- G2[, 5]; mAA <- G2[, 4]
  detX <- m11 * mLL - m1L * m1L
  singX <- !is.finite(detX) | detX == 0
  c1 <- ( mLL * m1A - m1L * mAL) / detX
  c2 <- (-m1L * m1A + m11 * mAL) / detX
  ZHZ <- mAA - (m1A * c1 + mAL * c2)
  dNg <- matrix(0, length(gs), 3); dNg[match(et, gs), ] <- dN
  ZHdN <- dNg[, 2] - (c1 * dNg[, 1] + c2 * dNg[, 3])
  ZHZ[singX] <- 0; ZHdN[singX] <- 0
  gam <- sum(ZHdN) / sum(ZHZ * dt)
  rcX <- (abs(detX) / pmax(m11, mLL)^2)[!singX]
  list(gamma = gam, mrd = surv0 - surv1,
       n_sing = sum(sing), n_singX = sum(singX),
       min_rc = if (length(rc_m)) min(rc_m) else NA_real_,
       min_rcX = if (length(rcX)) min(rcX) else NA_real_,
       flag = any(sing) || any(singX) || (length(rc_m) && min(rc_m) < RC_FLAG) ||
         (length(rcX) && min(rcX) < RC_FLAG))
}

# cell: one row of make_scenarios(); seed: all_seeds[rep]
paired_rep <- function(cell, rep, seed, mode = c("full", "advance"), B = 200L, t.hor = 1:5,
                       canon_nfail = NA_integer_) {
  mode <- match.arg(mode)
  t_start <- proc.time()[3]
  set.seed(seed)
  d <- dat_sim(n = cell$n, n.visit = DGM_FIXED$n.visit, gamma.0 = cell$gamma.0,
               gamma.L = cell$gamma.L, alpha.0 = cell$alpha.0, alpha.A = DGM_FIXED$alpha.A,
               alpha.L = DGM_FIXED$alpha.L, alpha.U = DGM_FIXED$alpha.U)
  ds0 <- ste_prepare_fast(d)
  ti0 <- tie_info(ds0)
  point_j <- NULL
  if (mode == "full" || ti0["nties"] > 0) {
    # canonical prefix, verbatim (consumes RNG only if the ORIGINAL data carry tied event times)
    est <- tryCatch(ste_estimate(d, n.visit = DGM_FIXED$n.visit, t.hor = t.hor,
                                 weight = "stab", estimand = "both"), error = function(e) NULL)
    if (is.null(est)) return(list(status = "point_failed", cell = cell$scenario, rep = rep))
    an <- tryCatch(ci_analytic(attr(est, "dat.seq")), error = function(e) NULL)
    point_j <- tryCatch(ste_estimate(d, n.visit = DGM_FIXED$n.visit, t.hor = t.hor,
                                     weight = "stab", estimand = "both"), error = function(e) NULL)
    if (is.null(point_j)) return(list(status = "point_failed", cell = cell$scenario, rep = rep))
    point_j <- list(gamma = point_j$gamma, mrd = point_j$mrd)
  }
  pe <- fit_exact_simul_guarded(ds0, d, t.hor)
  point_e <- list(gamma = pe$gamma, mrd = pe$mrd)

  rs <- make_resampler_idx(d)
  nt <- length(t.hor)
  J <- E1 <- E2 <- matrix(NA_real_, B, 1L + nt, dimnames = list(NULL, c("gamma", paste0("mrd", t.hor))))
  info <- matrix(NA_real_, B, 9L, dimnames = list(NULL, c("nties", "dt", "prep_fail", "jit_fail",
                   "e1_fail", "e2_fail", "n_sing", "n_singX", "min_rc")))
  e2_ties_left <- 0L
  for (b in seq_len(B)) {
    dl <- rs()
    dsb <- tryCatch(ste_prepare_fast(dl, n.visit = DGM_FIXED$n.visit), error = function(e) NULL)
    if (is.null(dsb)) { info[b, "prep_fail"] <- 1; next }
    info[b, "prep_fail"] <- 0
    tb <- tie_info(dsb); info[b, c("nties", "dt")] <- tb
    if (mode == "full") {
      fj <- tryCatch(ste_fit(dsb, dl, t.hor = t.hor, weight_col = "ipw.stab.cum", estimand = "both"),
                     error = function(e) NULL)
      info[b, "jit_fail"] <- is.null(fj)
      if (!is.null(fj)) J[b, ] <- c(fj$gamma, fj$mrd)
    } else {
      # advance the RNG exactly as the two aalen() calls' read.surv would
      if (tb["nties"] > 0) { stats::runif(tb["nties"]); stats::runif(tb["nties"]) }
    }
    e1 <- tryCatch(fit_exact_simul_guarded(dsb, dl, t.hor), error = function(e) NULL)
    info[b, "e1_fail"] <- is.null(e1)
    if (!is.null(e1)) {
      E1[b, ] <- c(e1$gamma, e1$mrd)
      info[b, c("n_sing", "n_singX", "min_rc")] <- c(e1$n_sing, e1$n_singX, e1$min_rc)
    }
    if (mode == "full" || b <= n_e2_check(cell$n) || is.null(e1) || isTRUE(e1$flag)) {
      e2 <- tryCatch(fit_exact_collapse(dsb, dl, attr(dl, "orig_id"), t.hor = t.hor),
                     error = function(e) NULL)
      info[b, "e2_fail"] <- is.null(e2)
      if (!is.null(e2)) { E2[b, ] <- c(e2$gamma, e2$mrd); e2_ties_left <- e2_ties_left + e2$ties_left }
    }
  }
  n_prep_fail <- sum(info[, "prep_fail"] == 1)
  aligned <- if (mode == "full") NA else (is.na(canon_nfail) || n_prep_fail == canon_nfail)
  list(status = "ok", cell = cell$scenario, scen_id = cell$scen_id, n = cell$n,
       alpha.0 = cell$alpha.0, gamma.0 = cell$gamma.0, gamma.L = cell$gamma.L,
       rep = rep, seed = seed, mode = mode, B = B, t.hor = t.hor,
       orig_nties = unname(ti0["nties"]),
       point_jitter = point_j, point_exact = point_e,
       J = J, E1 = E1, E2 = E2, info = info, e2_ties_left = e2_ties_left,
       n_prep_fail = n_prep_fail, canon_nfail = canon_nfail, advance_aligned = aligned,
       seconds = unname(proc.time()[3] - t_start),
       versions = c(R = as.character(getRversion()),
                    timereg = as.character(utils::packageVersion("timereg"))))
}
