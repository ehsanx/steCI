#' Simulate longitudinal survival data from an additive-hazard MSM
#'
#' Ports the data-generating mechanism of Keogh et al. (2021, Biometrical
#' Journal; 2023, Statistics in Medicine), as implemented in their public
#' `sequential_trials/dat_sim.R`. An unmeasured frailty `U` affects the
#' time-varying confounder `L` and the outcome hazard but not treatment `A`
#' (so it is not a confounder); treatment is monotone (once treated, always
#' treated). Event times come from a piecewise-constant generating hazard
#' `h`, with effective event hazard `max(h, 0)`: nonpositive values produce
#' no event in that interval.
#'
#' @param n Number of individuals.
#' @param n.visit Number of visits, K+1 (default 5): visits 0..K, admin
#'   censoring at t = n.visit.
#' @param gamma.0,gamma.L Treatment model `logit Pr(A=1) = gamma.0 + gamma.L*L`.
#'   `gamma.0` controls treatment prevalence; `gamma.L` controls confounding
#'   strength.
#' @param alpha.0,alpha.A,alpha.L,alpha.U Additive hazard
#'   `h = alpha.0 + alpha.A*A + alpha.L*L + alpha.U*U`. `alpha.0` controls the
#'   event rate; `alpha.A` is the conditional treatment coefficient in the
#'   generating hazard (default -0.04), distinct from the marginal E1 target.
#' @param seed Optional integer seed.
#'
#' @return A data.frame in long (start/stop) format with one row per
#'   at-risk visit interval and columns: `id`, `time` (interval start,
#'   0..K-1), `time.stop`, `event`, `visit`, `A`, `L`, `Alag1`, `L.baseline`.
#'   Also carries attribute `"neg_haz_frac"`: the fraction of hazard
#'   evaluations that were non-positive (relevant in low-event-rate regimes).
#' @references Keogh RH, Seaman SR, Gran JM, Vansteelandt S (2021).
#'   Simulating longitudinal data from marginal structural models using the
#'   additive hazard model. Biometrical Journal 63(7):1526-1541.
#' @export
dat_sim <- function(n = 1000, n.visit = 5,
                    gamma.0 = -1, gamma.L = 0.5,
                    alpha.0 = 0.2, alpha.A = -0.04,
                    alpha.L = 0.015, alpha.U = 0.015,
                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  expit <- function(x) exp(x) / (1 + exp(x))

  A <- matrix(nrow = n, ncol = n.visit)
  L <- matrix(nrow = n, ncol = n.visit)

  U <- stats::rnorm(n, 0, 0.1)
  L[, 1] <- stats::rnorm(n, U, 1)
  A[, 1] <- stats::rbinom(n, 1, expit(gamma.0 + gamma.L * L[, 1]))
  for (k in 2:n.visit) {
    L[, k] <- stats::rnorm(n, 0.8 * L[, k - 1] - A[, k - 1] + 0.1 * (k - 1) + U, 1)
    A[, k] <- ifelse(A[, k - 1] == 1, 1,
                     stats::rbinom(n, 1, expit(gamma.0 + gamma.L * L[, k])))
  }

  # Event times via inversion of a piecewise-constant additive hazard.
  T.obs <- rep(NA_real_, n)
  n_haz <- 0L; n_neg <- 0L
  for (k in 1:n.visit) {
    u.t <- stats::runif(n, 0, 1)
    haz <- alpha.0 + alpha.A * A[, k] + alpha.L * L[, k] + alpha.U * U
    n_haz <- n_haz + n
    n_neg <- n_neg + sum(haz <= 0, na.rm = TRUE)
    new.t <- -log(u.t) / haz
    # haz > 0 guard: additive hazard can go negative (low-event regimes).
    T.obs <- ifelse(is.na(T.obs) & new.t < 1 & haz > 0, k - 1 + new.t, T.obs)
  }
  D.obs <- ifelse(is.na(T.obs), 0, 1)
  T.obs <- ifelse(is.na(T.obs), n.visit, T.obs)

  visit_n <- n.visit - 1
  colnames(A) <- paste0("A.", 0:visit_n)
  colnames(L) <- paste0("L.", 0:visit_n)
  dat <- data.frame(id = 1:n, T.obs, D.obs, A, L)

  dat.long <- stats::reshape(
    data = dat,
    varying = c(paste0("A.", 0:visit_n), paste0("L.", 0:visit_n)),
    direction = "long", idvar = "id")
  dat.long <- dat.long[order(dat.long$id, dat.long$time), ]

  dat.long$time.stop <- dat.long$time + 1
  dat.long <- dat.long[dat.long$time < dat.long$T.obs, ]
  dat.long$time.stop <- ifelse(dat.long$time.stop > dat.long$T.obs,
                               dat.long$T.obs, dat.long$time.stop)
  dat.long$event <- ifelse(dat.long$time.stop == dat.long$T.obs & dat.long$D.obs == 1, 1, 0)
  dat.long$visit <- stats::ave(rep(1, nrow(dat.long)), dat.long$id, FUN = cumsum)

  dat.long <- dat.long %>%
    dplyr::group_by(.data$id) %>%
    dplyr::mutate(Alag1 = dplyr::lag(.data$A, n = 1),
                  Alag1 = tidyr::replace_na(.data$Alag1, 0),
                  L.baseline = dplyr::first(.data$L)) %>%
    dplyr::ungroup()

  out <- as.data.frame(dat.long)
  out <- out[, c("id", "time", "time.stop", "event", "visit",
                 "A", "L", "Alag1", "L.baseline")]
  attr(out, "neg_haz_frac") <- if (n_haz > 0) n_neg / n_haz else NA_real_
  out
}
