#' Delete-one-subject jackknife CIs for STE (M5)
#'
#' Leave-one-subject-out refitting of the full pipeline. Returns Wald CIs
#' from the jackknife standard error
#' `se = sqrt(((n-1)/n) * sum_k (theta_(k) - theta_hat)^2)` (centred at the
#' full-sample estimate, following Limozin et al.). Cost is `O(n)` refits, so
#' by default it is only run when the number of subjects is at most
#' `max_n`; this mirrors Limozin's restriction to small samples.
#'
#' @param dat.long Long data (subject-level `id`s).
#' @param n.visit,truncate,t.hor,weight,estimand See [ci_bootstrap()].
#' @param level Confidence level (default 0.95).
#' @param max_n Skip (with a warning) if the number of subjects exceeds this
#'   (default 500); set `Inf` to force.
#' @param ties,time.grid,point See [ci_bootstrap()]: the tie convention and
#'   grid for the full-sample estimate and every refit, and an optional
#'   full-sample estimate to centre on.
#'
#' @return A list with `point`, `gamma` and/or `mrd` CI data.frames, and
#'   `n_fail`.
#' @export
ci_jackknife <- function(dat.long, n.visit = 5, t.hor = 1:5,
                         weight = c("stab", "trunc"),
                         estimand = c("both", "hd", "mrd"),
                         level = 0.95, truncate = 0.95, max_n = 500,
                         ties = c("timereg", "simultaneous"), time.grid = NULL,
                         point = NULL) {
  weight <- match.arg(weight)
  estimand <- match.arg(estimand)
  ties <- match.arg(ties)
  z <- stats::qnorm(1 - (1 - level) / 2)

  if (is.null(point)) {
    point <- ste_estimate(dat.long, n.visit = n.visit, t.hor = t.hor,
                          weight = weight, estimand = estimand, truncate = truncate,
                          ties = ties, time.grid = time.grid)
  } else {
    .check_point(point, ties, estimand, t.hor)
  }
  ids <- unique(dat.long$id)
  n <- length(ids)
  if (n > max_n) {
    warning(sprintf("jackknife skipped: %d subjects > max_n = %d", n, max_n))
    return(list(point = point, gamma = NULL, mrd = NULL, n_fail = NA_integer_, skipped = TRUE))
  }
  jk <- .make_jackknifer(dat.long)

  gamma_jk <- rep(NA_real_, n)
  mrd_jk <- matrix(NA_real_, nrow = n, ncol = length(t.hor))
  n_fail <- 0L

  for (k in seq_len(n)) {
    dl <- jk$drop_k(k)
    est <- tryCatch(
      ste_estimate(dl, n.visit = n.visit, t.hor = t.hor,
                   weight = weight, estimand = estimand, truncate = truncate,
                   ties = ties, time.grid = time.grid),
      error = function(e) NULL)
    if (is.null(est)) { n_fail <- n_fail + 1L; next }
    if (!is.null(est$gamma)) gamma_jk[k] <- est$gamma
    if (!is.null(est$mrd)) mrd_jk[k, ] <- est$mrd
  }

  # The delete-one factor is (n-1)/n for a jackknife over n subjects. When some refits fail we
  # keep only the m finite pseudo-values, and rescaling by (m-1)/m would treat the survivors as
  # if they were the whole sample -- understating the SE whenever m < n. Use the design size n
  # and inflate by n/m to account for the missing deletions. The shipped runs have n_fail = 0,
  # so this branch has never been exercised and no published number changes; it is guarded
  # because a future caller with failures would otherwise get a silently optimistic SE.
  jk_se <- function(vals, center, n_design) {
    vals <- vals[is.finite(vals)]
    m <- length(vals)
    if (m <= 1L) return(NA_real_)
    sqrt(((n_design - 1) / n_design) * (n_design / m) * sum((vals - center)^2))
  }

  out <- list(point = point, n_fail = n_fail, skipped = FALSE, ties = ties)
  if (estimand %in% c("both", "hd")) {
    se <- jk_se(gamma_jk, point$gamma, length(ids))
    out$gamma <- data.frame(target = "hazard_diff", type = "jackknife",
                            estimate = point$gamma, se = se,
                            lower = point$gamma - z * se,
                            upper = point$gamma + z * se,
                            stringsAsFactors = FALSE)
  }
  if (estimand %in% c("both", "mrd")) {
    out$mrd <- do.call(rbind, lapply(seq_along(t.hor), function(j) {
      se <- jk_se(mrd_jk[, j], point$mrd[j], length(ids))
      data.frame(tau = t.hor[j], target = paste0("mrd_t", t.hor[j]),
                 type = "jackknife", estimate = point$mrd[j], se = se,
                 lower = point$mrd[j] - z * se, upper = point$mrd[j] + z * se,
                 stringsAsFactors = FALSE)
    }))
  }
  out
}
