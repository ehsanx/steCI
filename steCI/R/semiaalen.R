#' Constant-hazard-difference fit with simultaneous ties (internal)
#'
#' A from-scratch weighted fit of the semiparametric additive hazards model
#' `dN_i(t) = Y_i(t){X_i' dB(t) + Z_i' gamma dt} + dM_i(t)` with
#' `X = (1, L.baseline)` and `Z = A.baseline` (McKeague and Sasieni, 1994),
#' i.e. the model `timereg::aalen(Surv(time.new, time.stop.new, event) ~
#' const(A.baseline) + L.baseline, weights = w, robust = 1)` that [ste_fit()]
#' and [ci_analytic()] fit with `ties = "timereg"`:
#' \deqn{\hat\gamma = \Big\{\int Z^\top H Z\,dt\Big\}^{-1} \int Z^\top H\,dN,
#'   \qquad H = W - WX(X^\top W X)^{-1}X^\top W,}
#' with the `dt` integral taken over the grid of all distinct stop times.
#'
#' The only difference from timereg is the treatment of **tied event times**:
#' timereg breaks them by adding random jitter to every copy after the first
#' (`timereg:::read.surv`), so its fit consumes the random-number stream and
#' processes the copies sequentially; here all events at a time enter the same
#' increment, computed on the risk set just before that time. On data without
#' tied event times the estimate and all three variances equal timereg's to
#' machine precision (validated to ~1e-15 relative: `gamma`, the model-based
#' `var.gamma`, and `robvar.gamma` with row-level clusters and with subject
#' clusters, including timereg's `max.clust` quantile binning).
#'
#' Variances (all sandwich forms around \eqn{C = \int Z^\top H Z\,dt}):
#' model-based \eqn{C^{-1}\{\sum_t Z^\top H\,\mathrm{diag}(dN)\,H Z\}C^{-1}}
#' (timereg's `var.gamma`); robust \eqn{C^{-1}\{\sum_k e_k e_k^\top\}C^{-1}},
#' where \eqn{e_k} sums \eqn{\int (Z^\top H)_{\cdot r}\,d\hat M_r} over the rows
#' `r` of cluster `k` (rows by default -- timereg's `robvar.gamma` -- or
#' subjects when `clusters` is given).
#'
#' @param dat.seq Prepared person-trial data from [ste_prepare()].
#' @param weight_col Weight column name.
#' @param variance Compute the variances (default `TRUE`); `FALSE` returns the
#'   point estimate only (the refitting path of the resampling methods).
#' @param clusters Optional cluster labels (length `nrow(dat.seq)`) for the
#'   cluster-robust variance; `NULL` skips it.
#' @param max.clust timereg's `max.clust` (default 1000): with more clusters
#'   than this, clusters are merged into `max.clust` quantile bins exactly as
#'   timereg does. Applied to both the row-level and the cluster variance.
#' @return A list: `gamma`, `var_model`, `var_robust`, `var_cluster` (`NULL`
#'   when not computed), `intZHZ`, `intZHdN` (as in timereg), `n_singular`
#'   (grid intervals skipped -- no contribution -- because the risk-set Gram
#'   matrix \eqn{X^\top W X} could not be inverted by `solve()`, e.g. a tail
#'   risk set holding one subject). timereg does not document its treatment of
#'   such intervals and its output there depends on rounding, so agreement with
#'   timereg is claimed only where every risk set is invertible.
#' @references McKeague IW, Sasieni PD (1994) Biometrika 81(3):501-514.
#' @keywords internal
.semiaalen <- function(dat.seq, weight_col = "ipw.stab.cum", variance = TRUE,
                       clusters = NULL, max.clust = 1000) {
  w <- .get_w(dat.seq, weight_col)
  X <- cbind(1, dat.seq$L.baseline); Z <- cbind(dat.seq$A.baseline)
  p <- ncol(X); q <- ncol(Z); nr <- nrow(X)
  start <- dat.seq$time.new; stop <- dat.seq$time.stop.new
  ev <- dat.seq$event == 1
  # timereg's time grid for a model with const() terms: 0 and every distinct stop time.
  tg <- sort(unique(c(0, stop[stop > 0])))
  C <- matrix(0, q, q); ZHdN <- numeric(q); Vm <- matrix(0, q, q)
  n_singular <- 0L
  keep <- if (variance) vector("list", length(tg)) else NULL
  for (s in seq_along(tg)[-1L]) {
    t <- tg[s]; dt <- t - tg[s - 1L]
    R <- which(start < t & stop >= t)
    if (!length(R)) next
    XR <- X[R, , drop = FALSE]; ZR <- Z[R, , drop = FALSE]; wR <- w[R]
    Minv <- tryCatch(solve(crossprod(XR, wR * XR)), error = function(e) NULL)
    if (is.null(Minv)) { n_singular <- n_singular + 1L; next }
    ZtWX <- crossprod(ZR, wR * XR)
    M <- ZtWX %*% Minv                                   # q x p
    C <- C + dt * (crossprod(ZR, wR * ZR) - M %*% t(ZtWX))
    aR <- wR * (ZR - XR %*% t(M))                        # row r: (Z'H)_{., r}
    dNR <- as.numeric(stop[R] == t & ev[R])
    if (any(dNR > 0)) {
      ZHdN <- ZHdN + colSums(aR * dNR)
      if (variance) Vm <- Vm + crossprod(aR * dNR)
    }
    if (variance) keep[[s]] <- list(R = R, dt = dt, Minv = Minv, aR = aR, dNR = dNR)
  }
  Ci <- solve(C)
  gamma <- as.vector(Ci %*% ZHdN)
  out <- list(gamma = gamma, var_model = NULL, var_robust = NULL, var_cluster = NULL,
              intZHZ = C, intZHdN = ZHdN, n_singular = n_singular)
  if (!variance) return(out)

  # Row-level iid contributions e_r = sum_t (Z'H)_{., r} dM_r(t), with the fitted
  # dM_r = dN_r - X_r' dB - Z_r' gamma dt and dB = (X'WX)^{-1} X'W (dN - Z gamma dt).
  eps <- matrix(0, nr, q)
  for (s in seq_along(tg)[-1L]) {
    k <- keep[[s]]; if (is.null(k)) next
    XR <- X[k$R, , drop = FALSE]; ZR <- Z[k$R, , drop = FALSE]; wR <- w[k$R]
    r0 <- k$dNR - as.vector(ZR %*% gamma) * k$dt
    dB <- k$Minv %*% crossprod(XR, wR * r0)
    eps[k$R, ] <- eps[k$R, ] + k$aR * (r0 - as.vector(XR %*% dB))
  }
  sand <- function(cl) {
    E <- rowsum(eps, .timereg_clusters(cl, max.clust))
    Ci %*% crossprod(E) %*% t(Ci)
  }
  out$var_model <- Ci %*% Vm %*% t(Ci)
  out$var_robust <- sand(seq_len(nr) - 1L)             # timereg default: one cluster per row
  if (!is.null(clusters)) out$var_cluster <- sand(clusters)
  out
}

# timereg's cluster coding (read.surv + aalen): relabel to 0..(K-1) in sorted order and, when
# there are more than max.clust clusters, merge them into quantile bins (strict `<` test, so
# exactly max.clust clusters are NOT binned).
.timereg_clusters <- function(cl, max.clust = 1000) {
  K <- length(unique(cl))
  cl <- as.integer(factor(cl, labels = seq_len(K))) - 1L
  if (!is.null(max.clust) && max.clust < K) {
    qq <- unique(stats::quantile(cl, probs = seq(0, 1, by = 1 / max.clust)))
    cl <- as.integer(cut(cl, breaks = qq, include.lowest = TRUE)) - 1L
  }
  cl
}
