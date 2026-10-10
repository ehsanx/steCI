#' @keywords internal
#' @importFrom dplyr %>%
#' @importFrom rlang .data
NULL

# Fetch and validate the analysis weight column.
#
# `dat.seq[[weight_col]]` returns NULL for a name that is not a column of the
# data frame, and timereg::aalen() treats `weights = NULL` as its all-ones
# default. Without this guard an unrecognised `weight_col` therefore returned
# the UNWEIGHTED -- i.e. artificial-censoring-unadjusted, and so confounded --
# fit as an ordinary numeric result, with no error and no warning (measured at
# n = 400, seed 11: gamma = -0.07875 against the correct weighted -0.10301, a
# 24% attenuation). Every entry point that reads a weight column goes through
# this helper, so a mistyped name now fails loudly instead of silently
# changing the estimand. In .aalen_if() it also converts what would otherwise
# be a cryptic "non-conformable arguments" into a named error.
.get_w <- function(dat.seq, weight_col) {
  adm <- "the admissible values are \"ipw.stab.cum\" and \"ipw.stab.cum.trunc\"."
  if (!is.character(weight_col) || length(weight_col) != 1L || is.na(weight_col))
    stop("`weight_col` must be a single non-NA character column name; ", adm,
         call. = FALSE)
  if (!weight_col %in% names(dat.seq))
    stop("`weight_col` = \"", weight_col, "\" is not a column of `dat.seq`; ",
         adm, call. = FALSE)
  w <- dat.seq[[weight_col]]
  if (!is.numeric(w) || length(w) != nrow(dat.seq))
    stop("weight column \"", weight_col,
         "\" must be numeric and of length nrow(dat.seq).", call. = FALSE)
  if (any(!is.finite(w)) || any(w < 0))
    stop("weight column \"", weight_col,
         "\" must be non-missing, finite and non-negative.", call. = FALSE)
  w
}

# Validate `time.grid` (ste_prepare): a single positive finite number, and the long data's
# `time` / `time.stop` must lie on the 1/time.grid grid. The tolerance is in GRID UNITS
# (distinct grid points are 1 apart there), so it is keyed to the grid rather than to zero:
# a value is "on the grid" if it is within 1e-6 grid units of an integer.
.check_time_grid <- function(dat.long, time.grid) {
  if (!is.numeric(time.grid) || length(time.grid) != 1L || !is.finite(time.grid) ||
      time.grid <= 0)
    stop("`time.grid` must be a single positive finite number (grid units per visit).",
         call. = FALSE)
  off <- function(x) max(abs(x * time.grid - round(x * time.grid)))
  bad <- max(off(dat.long$time), off(dat.long$time.stop))
  if (bad > 1e-6)
    stop(sprintf(paste0("`time`/`time.stop` are not on a 1/%s grid (max offset %.3g grid ",
                        "units); use time.grid = NULL for continuous-time data."),
                 format(time.grid), bad), call. = FALSE)
  invisible(TRUE)
}

# Validate a caller-supplied full-sample estimate used to centre resampling intervals: it
# must come from the same tie convention (a point estimate from one estimator centring the
# refits of another is exactly the mixing this argument exists to prevent) and carry the
# requested estimands. An estimate without a `ties` element predates the switch and was
# necessarily timereg's.
.check_point <- function(point, ties, estimand, t.hor) {
  pt <- if (is.null(point$ties)) "timereg" else point$ties
  if (!identical(pt, ties))
    stop(sprintf("`point` was computed with ties = \"%s\" but ties = \"%s\" was requested.",
                 pt, ties), call. = FALSE)
  if (estimand %in% c("both", "hd") && !is.numeric(point$gamma))
    stop("`point` has no `gamma` for estimand = \"", estimand, "\".", call. = FALSE)
  if (estimand %in% c("both", "mrd") &&
      (!is.numeric(point$mrd) || length(point$mrd) != length(t.hor)))
    stop("`point$mrd` must have one value per horizon in `t.hor`.", call. = FALSE)
  invisible(TRUE)
}

# TRUE if a ties = "simultaneous" fit skipped any singular risk set (diagnostic count).
.any_singular <- function(est) {
  ns <- c(if (is.list(est$fit_hd)) est$fit_hd$n_singular,
          if (is.list(est$fit_mrd)) est$fit_mrd$n_singular)
  length(ns) > 0L && any(ns > 0L)
}

# Summarize bootstrap draws into percentile, basic/pivot, and normal CIs.
.boot_ci <- function(draws, point, level = 0.95, label = NA_character_) {
  a <- (1 - level) / 2
  draws <- draws[is.finite(draws)]
  se <- stats::sd(draws)
  q <- stats::quantile(draws, c(a, 1 - a), na.rm = TRUE, names = FALSE)
  z <- stats::qnorm(1 - a)
  data.frame(
    target = label,
    type   = c("percentile", "basic", "normal"),
    estimate = point,
    se = se,
    lower = c(q[1], 2 * point - q[2], point - z * se),
    upper = c(q[2], 2 * point - q[1], point + z * se),
    n_valid = length(draws),
    stringsAsFactors = FALSE)
}

# Build a fast subject-level resampler. Precomputes per-subject row indices
# once; each call resamples whole subjects with replacement and returns the
# long data with fresh contiguous ids. Index-based (subset once + rep ids),
# so it is O(total rows) rather than O(n^2) from rbind-ing many frames --
# essential at large n. RNG-identical to the old rbind version when ids = 1:n.
.make_resampler <- function(dat.long) {
  ids <- unique(dat.long$id)
  n <- length(ids)
  idx <- split(seq_len(nrow(dat.long)), factor(dat.long$id, levels = ids))
  lens <- lengths(idx)
  function() {
    j <- sample.int(n, n, replace = TRUE)
    dl <- dat.long[unlist(idx[j], use.names = FALSE), , drop = FALSE]
    dl$id <- rep.int(seq_len(n), lens[j])
    dl
  }
}

# Leave-one-subject-out long data (index-based).
.make_jackknifer <- function(dat.long) {
  ids <- unique(dat.long$id)
  n <- length(ids)
  idx <- split(seq_len(nrow(dat.long)), factor(dat.long$id, levels = ids))
  lens <- lengths(idx)
  list(n = n, drop_k = function(k) {
    keep <- seq_len(n)[-k]
    dl <- dat.long[unlist(idx[keep], use.names = FALSE), , drop = FALSE]
    dl$id <- rep.int(seq_len(n - 1L), lens[keep])
    dl
  })
}
