#' Standardize weighted-Aalen cumulative coefficients to marginal survival
#'
#' Empirical standardization (g-computation) of the additive-hazards fit to
#' "always treat" vs "never treat" marginal survival curves, averaging over
#' the observed baseline-`L` distribution. Ports the standardization block of
#' Keogh et al.'s `analysis_seq.R`.
#'
#' @param fit A fitted [timereg::aalen()] object from the formula
#'   `Surv(time.new, time.stop.new, event) ~ A.baseline + L.baseline`
#'   (fully time-varying coefficients; no `const`).
#' @param t.hor Numeric vector of horizon times.
#' @param L0 Numeric vector of baseline-`L` values (the population to
#'   standardize over).
#'
#' @return A list with `surv1`, `surv0` (survival under always/never treat)
#'   and `survdiff = surv1 - surv0` at each `t.hor`.
#' @keywords internal
mrd_from_aalen <- function(fit, t.hor, L0) {
  cum <- fit$cum  # columns: time, (Intercept), A.baseline, L.baseline
  f_int <- stats::stepfun(cum[, 1], c(0, cum[, 2]))
  f_A   <- stats::stepfun(cum[, 1], c(0, cum[, 3]))
  f_L   <- stats::stepfun(cum[, 1], c(0, cum[, 4]))
  surv1 <- vapply(t.hor, function(h)
    mean(exp(-(f_int(h) + f_A(h) + f_L(h) * L0))), numeric(1))
  surv0 <- vapply(t.hor, function(h)
    mean(exp(-(f_int(h) + f_L(h) * L0))), numeric(1))
  list(surv1 = surv1, surv0 = surv0, survdiff = surv1 - surv0)
}

# The same standardisation from cumulative coefficients already evaluated at the
# horizons (rows of B = (B0, BA, BL) at each t.hor), e.g. from .aalen_if(); arithmetic
# in the same order as mrd_from_aalen().
.mrd_from_B <- function(B, L0) {
  surv1 <- vapply(seq_len(nrow(B)), function(s)
    mean(exp(-(B[s, 1] + B[s, 2] + B[s, 3] * L0))), numeric(1))
  surv0 <- vapply(seq_len(nrow(B)), function(s)
    mean(exp(-(B[s, 1] + B[s, 3] * L0))), numeric(1))
  list(surv1 = surv1, surv0 = surv0, survdiff = surv1 - surv0)
}

#' Fit the sequential-trials additive-hazards estimator (point estimates)
#'
#' Lightweight point estimator used both directly and inside the resampling
#' CIs. Fits weighted [timereg::aalen()] models on the prepared person-trial
#' data and returns the causal estimands.
#'
#' @param dat.seq Output of [ste_prepare()].
#' @param dat.long The corresponding long data (for the baseline-`L`
#'   standardization population).
#' @param t.hor Horizon times for the marginal risk difference (default 1:5).
#' @param weight_col Name of the weight column (default "ipw.stab.cum";
#'   use "ipw.stab.cum.trunc" for truncated weights).
#' @param estimand One of "both", "hd" (constant hazard difference), "mrd"
#'   (marginal risk difference curve).
#' @param ties How tied event times are handled. `"timereg"` (default; the
#'   behaviour of every published simulation number) fits with
#'   [timereg::aalen()], which breaks tied event times by adding random
#'   jitter to every copy after the first -- so a fit on tied data consumes
#'   the random-number stream and processes the copies sequentially.
#'   `"simultaneous"` fits the same two models with tied events entering one
#'   increment on a common risk set: the constant hazard difference through
#'   [.semiaalen()] and the cumulative coefficients for the risk difference
#'   through [.aalen_if()] (the M6 fit). The two settings agree to machine
#'   precision on data without tied event times; `"simultaneous"` is
#'   deterministic on tied data. Resampling (subject bootstrap) always
#'   creates ties, because duplicated subjects share event times.
#'
#' @return A list. For "hd": `gamma` (constant hazard difference of treatment)
#'   and the fitted object `fit_hd`. For "mrd": `surv1`, `surv0`,
#'   `survdiff = surv1-surv0`, `mrd = surv0-surv1` (risk difference), `t.hor`,
#'   and `fit_mrd`. Also `ties`. Under `"simultaneous"` the fit objects are the
#'   lists returned by [.semiaalen()] and [.aalen_if()], whose `n_singular`
#'   counts risk sets skipped as singular.
#' @export
ste_fit <- function(dat.seq, dat.long, t.hor = 1:5,
                    weight_col = "ipw.stab.cum",
                    estimand = c("both", "hd", "mrd"),
                    ties = c("timereg", "simultaneous")) {
  estimand <- match.arg(estimand)
  ties <- match.arg(ties)
  w <- dat.seq[[weight_col]]
  out <- list(t.hor = t.hor, estimand = estimand, ties = ties)

  if (estimand %in% c("both", "hd")) {
    if (ties == "timereg") {
      fit_hd <- timereg::aalen(
        survival::Surv(time.new, time.stop.new, event) ~
          const(A.baseline) + L.baseline,
        data = dat.seq, n.sim = 0, weights = w)
      out$gamma <- unname(fit_hd$gamma[1, 1])
    } else {
      fit_hd <- .semiaalen(dat.seq, weight_col = weight_col, variance = FALSE)
      out$gamma <- fit_hd$gamma[1]
    }
    out$fit_hd <- fit_hd
  }

  if (estimand %in% c("both", "mrd")) {
    L0 <- dat.long$L[dat.long$visit == 1]
    if (ties == "timereg") {
      fit_mrd <- timereg::aalen(
        survival::Surv(time.new, time.stop.new, event) ~ A.baseline + L.baseline,
        data = dat.seq, n.sim = 0, weights = w)
      st <- mrd_from_aalen(fit_mrd, t.hor, L0)
    } else {
      fit_mrd <- .aalen_if(dat.seq, t.hor, weight_col = weight_col,
                           influence = FALSE, singular = "skip")
      st <- .mrd_from_B(fit_mrd$B, L0)
    }
    out$surv1 <- st$surv1
    out$surv0 <- st$surv0
    out$survdiff <- st$survdiff        # survival difference (always - never)
    out$mrd <- st$surv0 - st$surv1     # risk difference Pr(T<=t|a=1)-Pr(T<=t|a=0)
    out$fit_mrd <- fit_mrd
  }
  out
}

#' One-call convenience: prepare + fit from long data
#'
#' @inheritParams ste_prepare
#' @param t.hor,weight,estimand,ties See [ste_fit()]. `weight` is "stab" or "trunc".
#' @return The [ste_fit()] list, with the prepared data attached as
#'   attribute `"dat.seq"`.
#' @export
ste_estimate <- function(dat.long, n.visit = 5, t.hor = 1:5,
                         weight = c("stab", "trunc"),
                         estimand = c("both", "hd", "mrd"),
                         truncate = 0.95,
                         ties = c("timereg", "simultaneous"),
                         time.grid = NULL) {
  weight <- match.arg(weight)
  ties <- match.arg(ties)
  wcol <- if (weight == "stab") "ipw.stab.cum" else "ipw.stab.cum.trunc"
  ds <- ste_prepare(dat.long, n.visit = n.visit, truncate = truncate, time.grid = time.grid)
  res <- ste_fit(ds, dat.long, t.hor = t.hor, weight_col = wcol, estimand = estimand,
                 ties = ties)
  attr(res, "dat.seq") <- ds
  res
}
