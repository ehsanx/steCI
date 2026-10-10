#' Weighted Aalen additive fit with influence functions (internal)
#'
#' A from-scratch weighted least-squares fit of the Aalen additive hazards
#' model `dN_i(t) = Y_i(t) X_i(t)' dB(t) + dM_i(t)` with `X = (1, A.baseline,
#' L.baseline)`, returning the cumulative regression functions at the horizon
#' times together with per-row and per-subject (cluster) influence functions.
#' These influence functions drive the analytic delta-method standard errors
#' and the multiplier (wild) bootstrap in [ci_lef()].
#'
#' Tied event times are processed **simultaneously**: every event at time
#' `t` enters the same increment `dB(t)`, computed on the risk set just before
#' `t`. On data without tied event times the cumulative coefficients match
#' [timereg::aalen()] to machine precision (timereg instead breaks tied event
#' times by random jitter, so on tied data the two differ by design), and the
#' influence-function sandwich reproduces `robvar.cum` to within a few percent.
#'
#' @param dat.seq Prepared person-trial data from [ste_prepare()].
#' @param t.hor Increasing horizon times.
#' @param weight_col Weight column name.
#' @param dwdtheta Optional weight derivative (Approach 2, see [ci_lef()]).
#' @param influence If `FALSE`, return only the cumulative coefficients (the
#'   point-estimation path used by [ste_fit()] with `ties = "simultaneous"`);
#'   the arithmetic for `B` is identical either way.
#' @param singular What to do at an event time whose risk-set Gram matrix
#'   cannot be inverted (e.g. a late risk set with no `A.baseline` variation):
#'   `"error"` (default, the behaviour [ci_lef()] relies on) or `"skip"`, which
#'   leaves that increment at zero -- what [timereg::aalen()] was observed to do
#'   for an exactly singular risk set -- and counts it in `n_singular`.
#' @return A list: `B` (`length(t.hor) x 3` cumulative coefficients),
#'   `IFrow` (`n_row x 3 x n_t` row influence functions), `IFc`
#'   (`n_subject x 3 x n_t` cluster influence functions), `clusters`, `t.hor`,
#'   `n_singular`. With `influence = FALSE` the influence-function elements are
#'   `NULL`.
#' @keywords internal
.aalen_if <- function(dat.seq, t.hor, weight_col = "ipw.stab.cum", dwdtheta = NULL,
                      influence = TRUE, singular = c("error", "skip")) {
  singular <- match.arg(singular)
  X <- cbind(1, dat.seq$A.baseline, dat.seq$L.baseline)
  p <- ncol(X); nr <- nrow(X)
  w <- dat.seq[[weight_col]]; id <- dat.seq$id
  start <- dat.seq$time.new; stop <- dat.seq$time.stop.new; ev <- dat.seq$event == 1
  etimes <- sort(unique(stop[ev])); m <- length(etimes)
  nt <- length(t.hor)
  if (anyNA(t.hor) || is.unsorted(t.hor))
    stop("'t.hor' must be non-missing and non-decreasing", call. = FALSE)
  if (!influence && !is.null(dwdtheta))
    stop("`dwdtheta` requires influence = TRUE", call. = FALSE)
  do_A <- !is.null(dwdtheta); q <- if (do_A) ncol(dwdtheta) else 0L

  Bsnap <- matrix(0, nt, p)
  IFrow_snap <- if (influence) array(0, c(nr, p, nt)) else NULL
  IF <- if (influence) matrix(0, nr, p) else NULL
  Bacc <- numeric(p)
  Asnap <- if (do_A) array(0, c(p, q, nt)) else NULL
  Acc <- if (do_A) matrix(0, p, q) else NULL
  snap <- 1L; n_singular <- 0L

  snapshot <- function(s) {
    Bsnap[snap, ] <<- Bacc
    if (influence) IFrow_snap[, , snap] <<- IF
    if (do_A) Asnap[, , snap] <<- Acc
  }
  for (j in seq_len(m)) {
    tj <- etimes[j]
    # snapshot any horizon strictly before this event (value = accumulation so far)
    while (snap <= nt && t.hor[snap] < tj) { snapshot(snap); snap <- snap + 1L }
    inR <- start < tj & stop >= tj
    Xr <- X[inR, , drop = FALSE]; wr <- w[inR]
    if (singular == "error") {
      Dinv <- solve(crossprod(Xr, wr * Xr))
    } else {
      # solve() refuses exactly singular systems and those with reciprocal condition number
      # below .Machine$double.eps -- the same set on which singular = "error" stops.
      Dinv <- tryCatch(solve(crossprod(Xr, wr * Xr)), error = function(e) NULL)
      if (is.null(Dinv)) { n_singular <- n_singular + 1L; next }
    }
    isev <- inR & stop == tj & ev
    dNsum <- if (any(isev)) colSums(w[isev] * X[isev, , drop = FALSE]) else numeric(p)
    dB <- as.vector(Dinv %*% dNsum); Bacc <- Bacc + dB
    if (influence) {
      pred <- as.vector(Xr %*% dB)
      dMr  <- as.numeric(stop[inR] == tj & ev[inR]) - pred
      IF[inR, ] <- IF[inR, ] + t((Dinv %*% t(Xr)) * rep(wr * dMr, each = p))
      if (do_A) Acc <- Acc + Dinv %*% (t(Xr) %*% (dwdtheta[inR, , drop = FALSE] * dMr))
    }
  }
  # remaining horizons at/after the last event
  while (snap <= nt) { snapshot(snap); snap <- snap + 1L }

  if (!influence)
    return(list(B = Bsnap, IFrow = NULL, IFc = NULL, A = NULL, clusters = NULL,
                t.hor = t.hor, n_singular = n_singular))
  cl <- factor(id); K <- nlevels(cl)
  IFc <- array(0, c(K, p, nt))
  for (s in seq_len(nt)) IFc[, , s] <- rowsum(IFrow_snap[, , s], cl)
  list(B = Bsnap, IFrow = IFrow_snap, IFc = IFc, A = Asnap,
       clusters = levels(cl), t.hor = t.hor, n_singular = n_singular)
}
