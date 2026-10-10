# Tie-exact (simultaneous) refits of ste_fit()'s two weighted-Aalen fits, for the M4 tie experiment.
#
# WHY. timereg 2.0.7's aalen() calls timereg:::read.surv(), which breaks tied EVENT times with
#   time2[index] <- time2[index] + runif(nties) * dt,   dt = 0.5 * min positive |diff(time2)|
# (every copy after the first of a tied group is nudged UPWARD and then processed SEQUENTIALLY, the
# earlier copies having already left the risk set). A nonparametric bootstrap resample duplicates
# subjects, so every M4 refit on a resample carries ~one tied group per duplicated event row.
#
# Two independent tie-exact implementations are provided; both process a tied group SIMULTANEOUSLY
# (one risk set, summed dN), i.e. they give the estimator that the bootstrap's multinomial frequency
# weights define.
#
#  (E1) fit_exact_simul(): a from-scratch weighted least-squares Aalen fit (the .aalen_if()
#       recursion, without the influence functions), evaluated with SUFFIX SUMS of the 3x3 weighted
#       Gram matrix grouped by entry time, so no risk-set sum is ever formed by subtraction.
#         * MRD fit  (formula ~ A.baseline + L.baseline, all time-varying): dB(t) = M(t)^{-1} dN(t) at
#           each distinct event time t, M(t) = sum_{start<t<=stop} w x x', x = (1, A, L); B(tau) is
#           the sum over event times <= tau (right-continuous, as timereg's cum + stats::stepfun);
#           standardised over L0 exactly as mrd_from_aalen().
#         * HD fit (formula ~ const(A.baseline) + L.baseline): the McKeague-Sasieni estimator
#             gamma = [ sum_s ZHZ(g_s) (g_s - g_{s-1}) ]^{-1} sum_s Z'H(g_s) dN(g_s),
#           ZHZ = M_AA - M_AX M_XX^{-1} M_XA, Z'H dN = dN_A - M_AX M_XX^{-1} dN_X, X = (1, L), Z = A,
#           on timereg's semiparametric time grid g = {0} U unique(all stop times) (read.surv with
#           npar = FALSE uses ALL stop times, not only event times, so the dt-integral is exact for
#           piecewise-constant risk sets). Validated against timereg's intZHZ / intZHdN / gamma.
#  (E2) fit_exact_collapse(): collapse each bootstrap duplicate group to ONE row carrying weight
#       w x multiplicity and call timereg::aalen() itself. Because every quantity in both fits is a
#       sum over rows of w * f(row), this is algebraically the simultaneous fit, with timereg's own
#       numerics (including its handling of a singular risk-set Gram matrix). The collapsed data are
#       tie-free unless two DIFFERENT original subjects share an event time (checked and reported).
#
# Neither implementation consumes random numbers, so they can be interleaved with the canonical
# timereg refits without perturbing the canonical RNG stream.

# ---- resampler, RNG-identical to steCI:::.make_resampler(), but also returning the draw ----------
make_resampler_idx <- function(dat.long) {
  ids <- unique(dat.long$id)
  n <- length(ids)
  idx <- split(seq_len(nrow(dat.long)), factor(dat.long$id, levels = ids))
  lens <- lengths(idx)
  function() {
    j <- sample.int(n, n, replace = TRUE)
    dl <- dat.long[unlist(idx[j], use.names = FALSE), , drop = FALSE]
    dl$id <- rep.int(seq_len(n), lens[j])
    attr(dl, "orig_id") <- ids[j]          # new id k  <->  original subject ids[j[k]]
    dl
  }
}

# ---- suffix-sum risk-set sums --------------------------------------------------------------------
# For each query time t, sum over rows with start < t <= stop of the columns of P (a row-wise
# product matrix). Rows are grouped by their (few) distinct entry times; within a group, rows
# sorted by stop and REVERSE-cumulated, so a query reads the sum over {stop >= t} directly.
.risk_sums <- function(start, stop, P, qt) {
  out <- matrix(0, length(qt), ncol(P))
  for (s0 in sort(unique(start))) {
    g <- which(start == s0)
    o <- g[order(stop[g])]
    st <- stop[o]
    RS <- apply(P[o, , drop = FALSE], 2, function(v) rev(cumsum(rev(v))))
    if (!is.matrix(RS)) RS <- matrix(RS, nrow = length(o))
    RS <- rbind(RS, 0)                      # sentinel: nothing left
    pos <- findInterval(qt, st, left.open = TRUE) + 1L   # first index with stop >= t
    use <- qt > s0
    if (any(use)) out[use, ] <- out[use, ] + RS[pos[use], , drop = FALSE]
  }
  out
}

# 3x3 symmetric Gram entries, column order of P: 11,1A,1L,AA,AL,LL (x = (1, A, L)).
.gram_P <- function(w, A, L) cbind(w, w * A, w * L, w * A * A, w * A * L, w * L * L)

# Solve M z = b for many 3x3 symmetric M (rows of G, order 11,1A,1L,AA,AL,LL) via the adjugate.
# Returns list(z = k x 3, det = k, rc = crude reciprocal condition proxy).
.solve3 <- function(G, b) {
  a <- G[, 1]; bb <- G[, 2]; cc <- G[, 3]; d <- G[, 4]; e <- G[, 5]; f <- G[, 6]
  # M = [a b c; b d e; c e f]
  c11 <- d * f - e * e; c12 <- -(bb * f - cc * e); c13 <- bb * e - cc * d
  c22 <- a * f - cc * cc; c23 <- -(a * e - bb * cc); c33 <- a * d - bb * bb
  det <- a * c11 + bb * c12 + cc * c13
  z1 <- (c11 * b[, 1] + c12 * b[, 2] + c13 * b[, 3]) / det
  z2 <- (c12 * b[, 1] + c22 * b[, 2] + c23 * b[, 3]) / det
  z3 <- (c13 * b[, 1] + c23 * b[, 2] + c33 * b[, 3]) / det
  list(z = cbind(z1, z2, z3), det = det,
       rc = abs(det) / pmax(a, d, f)^3)
}

# ---- (E1) simultaneous, from scratch -------------------------------------------------------------
fit_exact_simul <- function(ds, dl, t.hor = 1:5, weight_col = "ipw.stab.cum",
                            estimand = c("both", "hd", "mrd"), refine = TRUE) {
  estimand <- match.arg(estimand)
  w <- ds[[weight_col]]
  A <- ds$A.baseline; L <- ds$L.baseline
  start <- ds$time.new; stop <- ds$time.stop.new; ev <- ds$event == 1
  P <- .gram_P(w, A, L)
  out <- list(t.hor = t.hor)

  # event sums by distinct event time: dN = sum_e w_e (1, A_e, L_e)
  et <- sort(unique(stop[ev]))
  fe <- factor(stop[ev], levels = et)
  dN <- rowsum(cbind(w[ev], w[ev] * A[ev], w[ev] * L[ev]), fe, reorder = TRUE)

  if (estimand %in% c("both", "mrd")) {
    G <- .risk_sums(start, stop, P, et)
    sol <- .solve3(G, dN)
    dB <- sol$z
    if (refine) {                      # one step of iterative refinement (cheap, 3x3)
      r <- dN - cbind(G[, 1] * dB[, 1] + G[, 2] * dB[, 2] + G[, 3] * dB[, 3],
                      G[, 2] * dB[, 1] + G[, 4] * dB[, 2] + G[, 5] * dB[, 3],
                      G[, 3] * dB[, 1] + G[, 5] * dB[, 2] + G[, 6] * dB[, 3])
      dB <- dB + .solve3(G, r)$z
    }
    Bc <- apply(dB, 2, cumsum)
    if (!is.matrix(Bc)) Bc <- matrix(Bc, nrow = 1)
    k <- findInterval(t.hor, et)       # number of event times <= tau
    Bh <- matrix(0, length(t.hor), 3)
    Bh[k > 0, ] <- Bc[k[k > 0], , drop = FALSE]
    L0 <- dl$L[dl$visit == 1]
    surv1 <- vapply(seq_along(t.hor), function(h)
      mean(exp(-(Bh[h, 1] + Bh[h, 2] + Bh[h, 3] * L0))), numeric(1))
    surv0 <- vapply(seq_along(t.hor), function(h)
      mean(exp(-(Bh[h, 1] + Bh[h, 3] * L0))), numeric(1))
    out$B <- Bh; out$surv1 <- surv1; out$surv0 <- surv0
    out$mrd <- surv0 - surv1
    out$min_rc_mrd <- min(sol$rc)
  }

  if (estimand %in% c("both", "hd")) {
    # timereg semiparametric grid: {0} U unique(all stop times) (read.surv, npar = FALSE)
    g <- sort(unique(c(0, stop[stop > 0])))
    gs <- g[-1]; dt <- diff(g)
    G <- .risk_sums(start, stop, P, gs)
    # partition with X = (1, L), Z = A: M_XX = [11 1L; 1L LL], M_XZ = (1A, AL), M_ZZ = AA
    m11 <- G[, 1]; m1L <- G[, 3]; mLL <- G[, 6]; m1A <- G[, 2]; mAL <- G[, 5]; mAA <- G[, 4]
    detX <- m11 * mLL - m1L * m1L
    c1 <- ( mLL * m1A - m1L * mAL) / detX      # c = M_XX^{-1} M_XZ
    c2 <- (-m1L * m1A + m11 * mAL) / detX
    ZHZ <- mAA - (m1A * c1 + mAL * c2)
    dNg <- matrix(0, length(gs), 3)
    mi <- match(et, gs); dNg[mi, ] <- dN
    ZHdN <- dNg[, 2] - (c1 * dNg[, 1] + c2 * dNg[, 3])
    intZHZ <- sum(ZHZ * dt); intZHdN <- sum(ZHdN)
    out$gamma <- intZHdN / intZHZ
    out$intZHZ <- intZHZ; out$intZHdN <- intZHdN
    out$min_rc_hd <- min(abs(detX) / pmax(m11, mLL)^2)
  }
  out
}

# ---- (E2) collapse duplicates to frequency weights, then timereg -----------------------------------
# orig_of_new: vector mapping new id (1..n) -> original subject id (attr "orig_id" of the resample).
fit_exact_collapse <- function(ds, dl, orig_of_new, t.hor = 1:5, weight_col = "ipw.stab.cum",
                               estimand = c("both", "hd", "mrd")) {
  estimand <- match.arg(estimand)
  oid <- orig_of_new[ds$id]
  key <- paste(oid, ds$trial, ds$time.new, sep = "|")
  first <- !duplicated(key)
  mult <- as.vector(table(factor(key, levels = unique(key))))   # in order of first appearance
  dc <- ds[first, , drop = FALSE]
  dc$wfreq <- dc[[weight_col]] * mult
  # sanity: the rows collapsed together must be identical in everything the fits read
  chk <- c("time.new", "time.stop.new", "event", "A.baseline", "L.baseline", weight_col)
  for (v in chk) {
    ref <- dc[[v]][match(key, key[first])]
    if (!isTRUE(all(ds[[v]] == ref))) stop("collapse: duplicates differ in ", v)
  }
  ties_left <- sum(duplicated(dc$time.stop.new[dc$event == 1]))
  out <- list(t.hor = t.hor, ties_left = ties_left, n_rows = nrow(dc))
  if (estimand %in% c("both", "hd")) {
    f <- timereg::aalen(survival::Surv(time.new, time.stop.new, event) ~
                          const(A.baseline) + L.baseline,
                        data = dc, n.sim = 0, weights = dc$wfreq)
    out$gamma <- unname(f$gamma[1, 1])
    out$intZHZ <- unname(f$intZHZ[1, 1]); out$intZHdN <- unname(f$intZHdN[1, 1])
  }
  if (estimand %in% c("both", "mrd")) {
    f <- timereg::aalen(survival::Surv(time.new, time.stop.new, event) ~ A.baseline + L.baseline,
                        data = dc, n.sim = 0, weights = dc$wfreq)
    L0 <- dl$L[dl$visit == 1]            # standardisation population WITH multiplicity
    st <- mrd_from_aalen(f, t.hor, L0)
    out$surv1 <- st$surv1; out$surv0 <- st$surv0; out$mrd <- st$surv0 - st$surv1
  }
  out
}

# number of tied event-time copies timereg would jitter, and its dt
tie_info <- function(ds) {
  s <- ds$time.stop.new; e <- ds$event == 1
  nt <- sum(duplicated(s[e]))
  dd <- abs(diff(s)); dt <- if (nt > 0) 0.5 * min(dd[dd > 0]) else NA_real_
  c(nties = nt, dt = dt)
}
