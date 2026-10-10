# A base-R re-implementation of steCI::ste_prepare() for the M4 tie experiment's compute budget.
#
# ste_prepare() spends ~90% of its time in dplyr grouped-mutate overhead (~0.7 s at n = 300,
# ~2 s at n = 1000 per call), and the M4 bootstrap calls it B = 200 times per replicate. This
# version performs the SAME arithmetic on the SAME rows in the SAME order, so every column the
# fits read (id, trial, time.new, time.stop.new, event, A.baseline, L.baseline, ipw.stab.cum,
# ipw.stab.cum.trunc) is required to be IDENTICAL (identical(), not all.equal) to ste_prepare()'s
# -- checked by validate_tie_exact.R on original data AND on bootstrap resamples -- and the
# timereg arm of the experiment is further checked to reproduce the canonical bootstrap
# output replicate by replicate. Row order matters beyond numerics here: timereg's tie-breaking
# jitter assigns its uniforms in row order and computes dt from consecutive-row differences.
ste_prepare_fast <- function(dat.long, n.visit = 5, truncate = 0.95) {
  dat.long <- dat.long[order(dat.long$id, dat.long$time), , drop = FALSE]
  blocks <- vector("list", n.visit)
  for (g in 0:(n.visit - 1L)) {
    d <- dat.long[dat.long$time >= g, , drop = FALSE]
    fi <- match(d$id, d$id)                         # first row of each id (d sorted by id, time)
    rownum <- seq_len(nrow(d)) - fi + 1L
    alb <- d$Alag1[fi]
    if (g > 0) {
      keep <- alb == 0
      d <- d[keep, , drop = FALSE]; rownum <- rownum[keep]
      fi <- match(d$id, d$id)
    }
    d$trial <- g
    d$rownum <- rownum
    d$time.new <- d$time - d$trial
    d$time.stop.new <- d$time.stop - d$trial
    d$Alag1.baseline <- d$Alag1[fi]
    d$A.baseline <- d$A[fi]
    d$L.baseline <- d$L[fi]
    blocks[[g + 1L]] <- d
  }
  ds <- do.call(rbind, blocks)
  rownames(ds) <- NULL
  N <- nrow(ds)
  same_next <- c(ds$id[-1] == ds$id[-N] & ds$trial[-1] == ds$trial[-N], FALSE)
  ds$A.lead1 <- ifelse(same_next, c(ds$A[-1], NA), NA_real_)
  ds$L.lead1 <- ifelse(same_next, c(ds$L[-1], NA), NA_real_)
  ds$Anext.equal.to.baseline <- ifelse(ds$A.lead1 == ds$A.baseline, 1, 0)

  dat.seq <- ds[ds$A == ds$A.baseline, , drop = FALSE]
  rownames(dat.seq) <- NULL

  ctrl <- dat.seq$A.baseline == 0
  predrows <- dat.seq$rownum != n.visit
  wt.mod.denom <- stats::glm(Anext.equal.to.baseline ~ L.lead1,
                             family = "binomial", data = dat.seq[ctrl, , drop = FALSE])
  dat.seq$wt.denom <- 1
  dat.seq$wt.denom[predrows] <- stats::predict(
    wt.mod.denom, newdata = dat.seq[predrows, , drop = FALSE], type = "response")
  dat.seq$wt.denom <- ifelse(dat.seq$A.baseline == 1, 1, dat.seq$wt.denom)
  wt.mod.num <- stats::glm(Anext.equal.to.baseline ~ as.factor(rownum) * L.baseline,
                           family = "binomial", data = dat.seq[ctrl, , drop = FALSE])
  dat.seq$wt.num <- 1
  dat.seq$wt.num[predrows] <- stats::predict(
    wt.mod.num, newdata = dat.seq[predrows, , drop = FALSE], type = "response")
  dat.seq$wt.num <- ifelse(dat.seq$A.baseline == 1, 1, dat.seq$wt.num)
  dat.seq$ipw.stab <- dat.seq$wt.num / dat.seq$wt.denom

  M <- nrow(dat.seq)
  first_in_grp <- c(TRUE, !(dat.seq$id[-1] == dat.seq$id[-M] & dat.seq$trial[-1] == dat.seq$trial[-M]))
  lagv <- c(NA_real_, dat.seq$ipw.stab[-M]); lagv[first_in_grp] <- NA_real_
  dat.seq$ipw.stab.lag <- ifelse(dat.seq$rownum == 1, 1, lagv)
  grp <- cumsum(first_in_grp)
  dat.seq$ipw.stab.cum <- unlist(lapply(split(dat.seq$ipw.stab.lag, grp), cumprod),
                                 use.names = FALSE)
  if (!is.na(truncate)) {
    pct <- stats::quantile(dat.seq$ipw.stab.cum, truncate, na.rm = TRUE)
    dat.seq$ipw.stab.cum.trunc <- pmin(dat.seq$ipw.stab.cum, pct)
  } else {
    dat.seq$ipw.stab.cum.trunc <- dat.seq$ipw.stab.cum
  }
  dat.seq
}

# Columns the fits (timereg and tie-exact) and the RNG-relevant tie logic read.
PREP_COLS <- c("id", "trial", "rownum", "time.new", "time.stop.new", "event", "A.baseline",
               "L.baseline", "ipw.stab.cum", "ipw.stab.cum.trunc", "L", "visit")

prep_identical <- function(a, b, cols = PREP_COLS) {
  if (nrow(a) != nrow(b)) return(FALSE)
  all(vapply(cols, function(v) identical(as.numeric(a[[v]]), as.numeric(b[[v]])), logical(1)))
}
