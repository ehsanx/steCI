#' Nonparametric (subject-level) bootstrap CIs for STE (M4)
#'
#' Resamples whole subjects with replacement and re-runs the entire pipeline
#' inside each replicate — trial construction, IPACW weight estimation, and
#' the weighted-Aalen fit — as recommended by Keogh et al. (2023) and
#' implemented (for pooled logistic) by Limozin et al. Returns percentile,
#' basic/pivot (`2*est - quantile`), and normal bootstrap CIs for the
#' constant hazard difference and for the marginal risk difference at each
#' horizon.
#'
#' Resampling duplicates subjects, and duplicated subjects share event times,
#' so every resample has tied event times: `ties` decides how the refits treat
#' them (see [ste_fit()]).
#'
#' @param dat.long Long data (subject-level `id`s).
#' @param n.visit,truncate,t.hor See [ste_prepare()] / [ste_fit()].
#' @param weight "stab" or "trunc".
#' @param estimand "both", "hd", or "mrd".
#' @param B Number of bootstrap resamples (default 500).
#' @param level Confidence level (default 0.95).
#' @param seed Optional integer seed.
#' @param ties,time.grid Passed to [ste_estimate()] for the full-sample
#'   estimate and every refit (defaults: `"timereg"`, `NULL`).
#' @param point Optional full-sample estimate (an [ste_estimate()] result for
#'   the same data and settings) to centre the intervals on, so that they are
#'   centred on exactly the estimate the caller reports; its `ties` must match.
#'   Default `NULL` computes it here (after `seed` is set).
#'
#' @return A list with `point` (the full-sample estimates), `gamma` (a
#'   data.frame of CIs for the hazard difference, if requested), `mrd` (a
#'   data.frame of CIs per horizon, if requested), `n_fail` (failed
#'   replicates), `B`, `ties`, and `n_singular` (under `"simultaneous"`, the
#'   number of refits in which at least one singular risk set was skipped;
#'   `NA` under `"timereg"`).
#' @export
ci_bootstrap <- function(dat.long, n.visit = 5, t.hor = 1:5,
                         weight = c("stab", "trunc"),
                         estimand = c("both", "hd", "mrd"),
                         B = 500, level = 0.95, truncate = 0.95, seed = NULL,
                         ties = c("timereg", "simultaneous"), time.grid = NULL,
                         point = NULL) {
  weight <- match.arg(weight)
  estimand <- match.arg(estimand)
  ties <- match.arg(ties)
  if (!is.null(seed)) set.seed(seed)

  if (is.null(point)) {
    point <- ste_estimate(dat.long, n.visit = n.visit, t.hor = t.hor,
                          weight = weight, estimand = estimand, truncate = truncate,
                          ties = ties, time.grid = time.grid)
  } else {
    .check_point(point, ties, estimand, t.hor)
  }

  resample <- .make_resampler(dat.long)

  gamma_boot <- rep(NA_real_, B)
  mrd_boot <- matrix(NA_real_, nrow = B, ncol = length(t.hor))
  n_fail <- 0L
  n_singular <- if (ties == "simultaneous") 0L else NA_integer_

  for (b in seq_len(B)) {
    dl <- resample()
    est <- tryCatch(
      ste_estimate(dl, n.visit = n.visit, t.hor = t.hor,
                   weight = weight, estimand = estimand, truncate = truncate,
                   ties = ties, time.grid = time.grid),
      error = function(e) NULL)
    if (is.null(est)) { n_fail <- n_fail + 1L; next }
    if (ties == "simultaneous" && .any_singular(est)) n_singular <- n_singular + 1L
    if (!is.null(est$gamma)) gamma_boot[b] <- est$gamma
    if (!is.null(est$mrd)) mrd_boot[b, ] <- est$mrd
  }

  out <- list(point = point, n_fail = n_fail, B = B, ties = ties, n_singular = n_singular)
  if (estimand %in% c("both", "hd")) {
    out$gamma <- .boot_ci(gamma_boot, point$gamma, level, label = "hazard_diff")
  }
  if (estimand %in% c("both", "mrd")) {
    out$mrd <- do.call(rbind, lapply(seq_along(t.hor), function(j)
      cbind(tau = t.hor[j],
            .boot_ci(mrd_boot[, j], point$mrd[j], level,
                     label = paste0("mrd_t", t.hor[j])))))
  }
  out
}
