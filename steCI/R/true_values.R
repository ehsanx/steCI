#' Approximate the true marginal survival curves and risk difference
#'
#' Simulates a large pseudo-RCT under the additive-hazard DGM: everyone
#' "always treated" (A = 1 at all visits) versus everyone "never treated"
#' (A = 0), sharing the same baseline frailty `U` and baseline confounder
#' `L0` across arms to reduce Monte-Carlo noise, then takes Kaplan-Meier
#' survival at the horizons. Ports Keogh et al.'s `true_values.R`.
#'
#' @param n Pseudo-RCT size per arm (default 1e6).
#' @param n.visit,t.hor,alpha.0,alpha.A,alpha.L,alpha.U DGM settings (must
#'   match the [dat_sim()] cell whose truth you want). Treatment-model
#'   parameters are irrelevant here (treatment is set by regime).
#' @param seed Optional integer seed.
#'
#' @param n_hd Sub-sample size per arm used to estimate the causal marginal
#'   constant hazard difference (`gamma_true`) via a constant-effect additive
#'   model; the full `n` is used for the Kaplan-Meier survival curves.
#' @return A list: `t.hor`, `surv1`, `surv0`, `survdiff = surv1-surv0`,
#'   `mrd = surv0-surv1` (risk difference), and `gamma_true` --- the causal
#'   marginal constant hazard difference of "always treat" vs "never treat".
#'   These are Monte-Carlo approximations to the single-cohort marginal
#'   targets ("always treated" versus "never treated" from baseline). The
#'   marginal hazard difference differs from the DGM's *conditional*
#'   coefficient `alpha.A` because sustained treatment also lowers the
#'   time-varying confounder L. The paper's simulation study scores against
#'   the closed-form values stored in the compendium's
#'   `dev/truth_cache_exact.rds` (built by `dev/truth_closed_form.R`).
#' @export
true_values <- function(n = 1e6, n.visit = 5, t.hor = 1:5,
                        alpha.0 = 0.2, alpha.A = -0.04,
                        alpha.L = 0.015, alpha.U = 0.015, seed = NULL,
                        n_hd = 2e5) {
  if (!is.null(seed)) set.seed(seed)
  U  <- stats::rnorm(n, 0, 0.1)
  L0 <- stats::rnorm(n, U, 1)

  arm <- function(aval) {
    L <- matrix(nrow = n, ncol = n.visit)
    A <- matrix(aval, nrow = n, ncol = n.visit)
    L[, 1] <- L0
    for (k in 2:n.visit) {
      L[, k] <- stats::rnorm(n, 0.8 * L[, k - 1] - A[, k - 1] + 0.1 * (k - 1) + U, 1)
    }
    Tt <- rep(NA_real_, n)
    for (k in 1:n.visit) {
      u <- stats::runif(n, 0, 1)
      haz <- alpha.0 + alpha.A * A[, k] + alpha.L * L[, k] + alpha.U * U
      nt <- -log(u) / haz
      Tt <- ifelse(is.na(Tt) & nt < 1 & haz > 0, k - 1 + nt, Tt)
    }
    D <- ifelse(is.na(Tt), 0, 1)
    Tt <- ifelse(is.na(Tt), n.visit, Tt)
    list(Tt = Tt, D = D)
  }

  a1 <- arm(1); a0 <- arm(0)
  km1 <- survival::survfit(survival::Surv(a1$Tt, a1$D) ~ 1)
  km0 <- survival::survfit(survival::Surv(a0$Tt, a0$D) ~ 1)
  s1 <- summary(km1, times = t.hor, extend = TRUE)$surv
  s0 <- summary(km0, times = t.hor, extend = TRUE)$surv

  # Causal marginal constant hazard difference: constant-effect additive model
  # on a randomized (always- vs never-treat) sub-sample.
  m <- min(n_hd, n); idx <- seq_len(m)
  dd <- data.frame(Tstop = c(a1$Tt[idx], a0$Tt[idx]),
                   event = c(a1$D[idx], a0$D[idx]),
                   A = rep(c(1, 0), each = m))
  gamma_true <- timereg::aalen(survival::Surv(Tstop, event) ~ const(A),
                               data = dd, n.sim = 0)$gamma[1, 1]

  list(t.hor = t.hor, surv1 = s1, surv0 = s0,
       survdiff = s1 - s0, mrd = s0 - s1, gamma_true = unname(gamma_true))
}
