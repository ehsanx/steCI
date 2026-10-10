# Shared by applied/heart_analysis.R (Tables heart, heart_mrd; Section 5 prose) and
# dev/heart_grid_sensitivity.R (Web Table S9). Everything a printed heart number depends on
# that the two scripts must agree on lives here, once:
#
#   * heart_to_ste()  -- the jasa -> person-period mapping (was duplicated in both scripts);
#   * HEART_TIES      -- THE tie-handling switch for every printed heart estimate;
#   * HEART_SEED      -- set before the first random draw of either script.
#
# ---- THE TIE SWITCH -------------------------------------------------------------------------
# timereg::aalen() breaks tied event times by adding random jitter to every copy after the
# first (timereg:::read.surv), so a fit on tied data consumes the random-number stream and
# processes the copies sequentially, not simultaneously. The heart data are on a 1-day grid
# and have many exact ties; every subject bootstrap resample adds more (duplicated subjects).
#   "simultaneous": tied events enter one increment on a common risk set (the M6/.aalen_if
#                   convention). Point estimates, M1-M3, the M4/M5 refits and their centring
#                   are all tie-exact and consume no random numbers; the MRD point estimate
#                   is then the M6 fit itself.
#   "timereg":      timereg's jittered fits everywhere, including the point estimates (the
#                   behaviour of the published simulation). M6 standard errors still come
#                   from the simultaneous-tie .aalen_if linearisation, which has no jitter.
# The two settings agree to machine precision on data without tied event times.
# Convention of the paper: "simultaneous" for every printed heart
# number; the primary simulation keeps the timereg convention. Override per run with the
# environment variable HEART_TIES=simultaneous|timereg.
HEART_TIES_DEFAULT <- "simultaneous"   # the paper's heart convention; the primary simulation keeps timereg
HEART_SEED <- 2024L

heart_ties <- function() {
  x <- Sys.getenv("HEART_TIES", unset = "")
  if (!nzchar(x)) x <- HEART_TIES_DEFAULT
  if (is.na(x))
    stop("No tie convention chosen: set HEART_TIES=simultaneous or HEART_TIES=timereg ",
         "(or HEART_TIES_DEFAULT in applied/heart_common.R).", call. = FALSE)
  match.arg(x, c("simultaneous", "timereg"))
}

# Seed with the RNG kinds pinned, so the result does not depend on whatever RNGkind() the
# calling session happened to have (the value of set.seed(s) alone depends on the kind).
heart_seed <- function() {
  set.seed(HEART_SEED, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
}

# ---- map jasa to the steCI id/time/time.stop/event/A/L schema ----
# Person-period expansion: w-day visits, horizon K visits. A = transplanted by the
# start of the visit (absorbing); L = standardized age (baseline confounder of both
# transplant priority and mortality); event = death within the visit. Times are in visit
# units on a 1/w grid (whole days), which is why every ste_prepare() call on these data
# passes time.grid = w: it makes events on the same trial-relative day exactly tied.
heart_to_ste <- function(w = 30, K = 6) {
  jasa <- survival::jasa
  n_all <- nrow(jasa)
  jasa <- jasa[jasa$futime > 0, ]                    # drop zero-follow-up (died day 0)
  if (nrow(jasa) < n_all)
    message(sprintf("dropped %d patient(s) with non-positive follow-up", n_all - nrow(jasa)))
  n <- nrow(jasa); Lage <- as.numeric(scale(jasa$age))
  out <- vector("list", n)
  for (i in seq_len(n)) {
    fut <- jasa$futime[i]; dead <- jasa$fustat[i] == 1; wt <- jasa$wait.time[i]
    rr <- list()
    for (k in 0:(K - 1)) {
      day0 <- k * w; day1 <- (k + 1) * w
      if (fut <= day0) break                         # left before this visit
      Ak <- as.integer(!is.na(wt) & wt <= day0)      # transplanted by visit start
      stop_day <- min(day1, fut)
      ev <- as.integer(dead & fut > day0 & fut <= day1)
      rr[[length(rr) + 1]] <- data.frame(
        id = i, time = k, time.stop = k + (stop_day - day0) / w,
        event = ev, A = Ak, L = Lage[i])
    }
    out[[i]] <- do.call(rbind, rr)
  }
  d <- do.call(rbind, out)
  d <- d[order(d$id, d$time), ]
  d$visit <- ave(rep(1, nrow(d)), d$id, FUN = cumsum)
  d <- d %>% dplyr::group_by(id) %>%
    dplyr::mutate(Alag1 = dplyr::lag(A, 1), Alag1 = tidyr::replace_na(Alag1, 0),
                  L.baseline = dplyr::first(L)) %>% dplyr::ungroup()
  as.data.frame(d)
}

# Event-time near-ties in prepared data: pairs of distinct event times closer than `tol`
# (visit units). Must be 0 once time.grid is used; exact ties are fine.
heart_near_ties <- function(dat.seq, tol = 1e-9) {
  et <- sort(unique(dat.seq$time.stop.new[dat.seq$event == 1]))
  sum(diff(et) < tol)
}
