#' Build sequentially emulated trials with stabilized IPACW weights
#'
#' Ports `analysis_seq.R` (Keogh et al. 2023, public `sequential_trials`
#' repo): from long data it constructs nested trials with origins
#' 0..(n.visit-1), restricts each later trial to individuals untreated at
#' the previous visit, and imposes artificial censoring when control-arm
#' subjects initiate treatment. This implementation requires absorbing
#' observed treatment (once treated, always treated). Its two stabilized
#' inverse-probability-of-artificial-censoring weight (IPACW) models are
#' fitted to the control arm; treated-arm weights are one.
#'
#' @param dat.long Long data from [dat_sim()] (or same schema): columns
#'   `id, time, time.stop, event, visit, A, L, Alag1, L.baseline`, with one
#'   start-stop record per subject-interval. Row order is part of the
#'   contract: every derived quantity below (`rownum`, `first()`, `lead()`,
#'   `lag()`) is positional in within-subject row order rather than keyed on
#'   `time`, so the function sorts its input by `(id, time)` on entry.
#' @param n.visit Number of visits (default 5).
#' @param truncate Upper quantile for weight truncation (default 0.95);
#'   set `NA` to skip.
#' @param time.grid Optional number of grid units per visit for data recorded
#'   on a discrete sub-visit grid (e.g. `30` for daily data on 30-day visits,
#'   so every `time.stop` is a multiple of `1/30`). The trial-relative stop
#'   time `time.stop - trial` is otherwise formed in floating point, and two
#'   events at the same trial-relative day reached from different trial
#'   origins can then differ in the last bit (e.g. `(1 + 11/30) - 1` versus
#'   `11/30`), so that tied event times are processed as distinct, in an
#'   arbitrary rounding-determined order. With `time.grid` set,
#'   `time.stop.new` is snapped to the exact grid value `round(x * time.grid)
#'   / time.grid`, which makes such ties exact (`time.new`, a difference of
#'   visit indices, is already exact); the input must lie on the grid
#'   (checked to `1e-6` grid units). Default `NULL` leaves the times exactly
#'   as before, which is what continuous-time data (e.g. [dat_sim()]) require.
#'
#' @return A data.frame of stacked, artificially-censored person-trial rows
#'   with the trial-relative time columns `time.new`, `time.stop.new`, the
#'   fixed `A.baseline`, `L.baseline`, and analysis weights `ipw.stab.cum`
#'   (and `ipw.stab.cum.trunc`). Ready for [ste_fit()].
#' @export
ste_prepare <- function(dat.long, n.visit = 5, truncate = 0.95, time.grid = NULL) {
  if (!is.null(time.grid)) .check_time_grid(dat.long, time.grid)
  # Enforce the ordering invariant the body assumes: rownum/first()/lead()/lag()
  # are positional within subject, not keyed on time. This is a no-op on
  # dat_sim() output (already sorted at dat_sim.R:75) and on bootstrap /
  # jackknife resamples (.make_resampler and .make_jackknifer in utils.R
  # concatenate whole per-subject index blocks, preserving within-subject
  # order), so no published number changes.
  dat.long <- dat.long[order(dat.long$id, dat.long$time), , drop = FALSE]

  origins <- 0:(n.visit - 1)

  build_trial <- function(g) {
    d <- dat.long[dat.long$time >= g, , drop = FALSE]
    d <- d %>%
      dplyr::mutate(trial = g) %>%
      dplyr::group_by(.data$id) %>%
      dplyr::mutate(rownum = dplyr::row_number(),
                    time.new = .data$time - .data$trial,
                    time.stop.new = .data$time.stop - .data$trial,
                    Alag1.baseline = dplyr::first(.data$Alag1))
    if (g > 0) {
      d <- dplyr::filter(d, .data$Alag1.baseline == 0)  # untreated at trial entry
    }
    d %>%
      dplyr::mutate(A.baseline = dplyr::first(.data$A),
                    L.baseline = dplyr::first(.data$L)) %>%
      dplyr::ungroup()
  }

  dat.seq.orig <- dplyr::bind_rows(lapply(origins, build_trial))
  if (!is.null(time.grid)) {
    # Exact grid values: every trial-relative stop time becomes fl(N / time.grid) for an
    # integer N, so equal grid times are bit-identical whichever trial origin they were
    # reached from. (time.new = time - trial is a difference of visit indices: already exact.)
    dat.seq.orig$time.stop.new <- round(dat.seq.orig$time.stop.new * time.grid) / time.grid
  }

  # Indicator that the *next* treatment still equals trial-baseline treatment.
  dat.seq.orig <- dat.seq.orig %>%
    dplyr::group_by(.data$id, .data$trial) %>%
    dplyr::mutate(A.lead1 = dplyr::lead(.data$A, n = 1),
                  L.lead1 = dplyr::lead(.data$L, n = 1)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Anext.equal.to.baseline = ifelse(.data$A.lead1 == .data$A.baseline, 1, 0))

  # Artificial censoring: keep rows where current A equals baseline A.
  dat.seq <- dat.seq.orig[dat.seq.orig$A == dat.seq.orig$A.baseline, , drop = FALSE]
  dat.seq <- as.data.frame(dat.seq)

  ctrl <- dat.seq$A.baseline == 0
  predrows <- dat.seq$rownum != n.visit  # last overall interval carries no next-step weight

  # Denominator model (conditions on next/current time-varying L).
  wt.mod.denom <- stats::glm(Anext.equal.to.baseline ~ L.lead1,
                             family = "binomial", data = dat.seq[ctrl, , drop = FALSE])
  dat.seq$wt.denom <- 1
  dat.seq$wt.denom[predrows] <- stats::predict(
    wt.mod.denom, newdata = dat.seq[predrows, , drop = FALSE], type = "response")
  dat.seq$wt.denom <- ifelse(dat.seq$A.baseline == 1, 1, dat.seq$wt.denom)

  # Numerator model (conditions on trial-baseline L and visit number).
  wt.mod.num <- stats::glm(Anext.equal.to.baseline ~ as.factor(rownum) * L.baseline,
                           family = "binomial", data = dat.seq[ctrl, , drop = FALSE])
  dat.seq$wt.num <- 1
  dat.seq$wt.num[predrows] <- stats::predict(
    wt.mod.num, newdata = dat.seq[predrows, , drop = FALSE], type = "response")
  dat.seq$wt.num <- ifelse(dat.seq$A.baseline == 1, 1, dat.seq$wt.num)

  dat.seq$ipw.stab <- dat.seq$wt.num / dat.seq$wt.denom

  # Analysis weight = cumulative product of LAGGED stabilized ratios (row 1 = 1).
  dat.seq <- dat.seq %>%
    dplyr::group_by(.data$id, .data$trial) %>%
    dplyr::mutate(ipw.stab.lag = dplyr::lag(.data$ipw.stab, n = 1),
                  ipw.stab.lag = ifelse(.data$rownum == 1, 1, .data$ipw.stab.lag),
                  ipw.stab.cum = cumprod(.data$ipw.stab.lag)) %>%
    dplyr::ungroup()
  dat.seq <- as.data.frame(dat.seq)

  if (!is.na(truncate)) {
    pct <- stats::quantile(dat.seq$ipw.stab.cum, truncate, na.rm = TRUE)
    dat.seq$ipw.stab.cum.trunc <- pmin(dat.seq$ipw.stab.cum, pct)
  } else {
    dat.seq$ipw.stab.cum.trunc <- dat.seq$ipw.stab.cum
  }
  dat.seq
}
