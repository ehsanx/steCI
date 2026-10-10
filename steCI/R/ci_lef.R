#' Analytic delta-method and multiplier/LEF bootstrap CIs for the MRD (M6)
#'
#' Confidence intervals for the marginal risk-difference curve based on the
#' influence functions of the weighted Aalen additive-hazards estimator
#' ([.aalen_if()]). This adapts the linearised-estimating-function (LEF)
#' bootstrap of Limozin, Seaman & Su (2025) to the weighted-Aalen
#' additive-hazards estimator: because that estimator is weighted least squares, its
#' influence-function linearisation is asymptotically (not exactly) linear --- the
#' Gram matrix is data-dependent --- and the bootstrap reduces to a fast **multiplier
#' (wild) bootstrap** of the influence functions --- no refitting.
#'
#' Provides, at each horizon: (i) analytic delta-method Wald CIs from the
#' row-level and **cluster-robust** sandwich standard errors of the MRD (the
#' cluster version accounts for a subject appearing in several stacked
#' trials); (ii) pointwise multiplier-bootstrap percentile CIs; and (iii) a
#' **simultaneous confidence band** for the whole MRD curve (sup-t).
#'
#' By default the multiplier bootstrap treats the estimated IPACW weights and the
#' standardisation population as fixed (an Approach-1 analogue); the
#' nonparametric bootstrap ([ci_bootstrap()]) relaxes both. Set
#' `standardization = "random"` to additionally propagate the sampling
#' variability of the empirical baseline-covariate (standardisation) population,
#' and `approach = "a2"` to propagate weight-estimation uncertainty.
#' This first-order correction was evaluated empirically; its asymptotic
#' validity is not established in the paper.
#'
#' @param dat.long Long data (subject-level `id`s).
#' @param n.visit,truncate,t.hor,weight See [ci_bootstrap()].
#' @param B Number of multiplier draws (default 500).
#' @param level Confidence level (default 0.95).
#' @param multiplier Multiplier distribution: "rademacher" (default),
#'   "normal", or "poisson" (mean-centred unit-variance).
#' @param approach Weight-uncertainty handling: `"a1"` (default) treats the
#'   estimated IPACW weights as fixed; `"a2"` additionally propagates
#'   weight-estimation uncertainty through the two logistic weight models.
#' @param standardization Standardisation-population handling for the marginal
#'   risk difference: `"fixed"` (default) conditions on the observed baseline
#'   covariates; `"random"` adds the empirical-process (g-computation) term
#'   \eqn{(S_0^{(i)}-S_1^{(i)})-\widehat{MRD}} to each subject's influence
#'   function, giving the *unconditional* variance for the population estimand.
#'   Affects only the subject-level methods (cluster analytic, multiplier, band).
#' @param seed Optional integer seed.
#' @param time.grid Passed to [ste_prepare()] (default `NULL`). The fit
#'   ([.aalen_if()]) always processes tied event times simultaneously; for
#'   data on a discrete sub-visit grid, `time.grid` makes those ties exact
#'   rather than last-bit near-ties.
#'
#' @return A list with `point` (MRD at each horizon), `mrd` (a data.frame of
#'   pointwise CIs: methods `analytic_row`, `analytic_cluster`, `multiplier`),
#'   and `band` (a data.frame with the simultaneous band and its critical
#'   value `calpha`).
#' @references Hu F, Kalbfleisch JD (2000) Canadian Journal of Statistics
#'   28(3):449-499. Lin DY, Ying Z (1994) Biometrika 81(1):61-71.
#' @export
ci_lef <- function(dat.long, n.visit = 5, t.hor = 1:5,
                   weight = c("stab", "trunc"), B = 500, level = 0.95,
                   multiplier = c("rademacher", "normal", "poisson"),
                   approach = c("a1", "a2"),
                   standardization = c("fixed", "random"),
                   truncate = 0.95, seed = NULL, time.grid = NULL) {
  weight <- match.arg(weight)
  multiplier <- match.arg(multiplier)
  approach <- match.arg(approach)
  standardization <- match.arg(standardization)
  if (!is.null(seed)) set.seed(seed)
  wcol <- if (weight == "stab") "ipw.stab.cum" else "ipw.stab.cum.trunc"
  # Approach 2 differentiates the UNTRUNCATED cumulative stabilised weight
  # w_r = prod_{s<r} p_num(s)/p_denom(s) (Web Appendix B, Eq. (S9)). Under
  # weight = "trunc" the analysis weight is ipw.stab.cum.trunc, whose derivative is zero
  # wherever truncation binds, so pairing it with the Eq. (S9) derivative would combine
  # one weight with another weight's derivative. Refuse rather than do that silently.
  if (approach == "a2" && weight == "trunc")
    stop("approach = \"a2\" is derived for the untruncated stabilised weight ",
         "(Web Appendix B, Eq. (S9)); use weight = \"stab\".")

  ds <- ste_prepare(dat.long, n.visit = n.visit, truncate = truncate, time.grid = time.grid)
  L0 <- dat.long$L[dat.long$visit == 1]

  # Approach 2 additionally propagates weight-estimation uncertainty: compute the
  # weight-derivative dw/dtheta and the weight-model influence functions, and the
  # sensitivity A(t) = dB(t)/dtheta, then add A(t) %*% IF_theta_c to the cluster
  # influence function.
  wif <- if (approach == "a2") .weight_if(ds, n.visit) else NULL
  fit <- .aalen_if(ds, t.hor, weight_col = wcol,
                   dwdtheta = if (approach == "a2") wif$dwdtheta else NULL)

  nt <- length(t.hor); K <- dim(fit$IFc)[1]
  z <- stats::qnorm(1 - (1 - level) / 2)
  a <- (1 - level) / 2
  # baseline L per cluster (subject), aligned to the rows of fit$IFc / fit$clusters
  Lc <- ds$L.baseline[match(fit$clusters, as.character(ds$id))]

  mrd <- numeric(nt); se_row <- numeric(nt); se_cl <- numeric(nt)
  IFc_mrd <- matrix(0, K, nt)
  for (s in seq_len(nt)) {
    B0 <- fit$B[s, 1]; BA <- fit$B[s, 2]; BL <- fit$B[s, 3]
    S1 <- exp(-(B0 + BA + BL * L0)); S0 <- exp(-(B0 + BL * L0))
    mrd[s] <- mean(S0 - S1)                               # risk difference
    g <- c(mean(S1 - S0), mean(S1), mean(L0 * (S1 - S0))) # gradient wrt (B0,BA,BL)
    IFc_B <- fit$IFc[, , s]
    if (approach == "a2") IFc_B <- IFc_B + t(fit$A[, , s] %*% t(wif$IFtheta))
    IFc_mrd[, s] <- IFc_B %*% g
    if (standardization == "random") {
      # empirical-process (g-computation) term for the random standardisation
      # population. Each subject contributes (S0_i - S1_i - MRD)/K: the 1/K matches
      # the unnormalised convention (fit$IFc rows are O(1/K)), since MRD is itself
      # the sample mean (1/K) sum_i (S0_i - S1_i). Centred to be mean-zero across
      # clusters so it is a proper influence-function addend.
      S1c <- exp(-(B0 + BA + BL * Lc)); S0c <- exp(-(B0 + BL * Lc))
      mc <- S0c - S1c
      IFc_mrd[, s] <- IFc_mrd[, s] + (mc - mean(mc)) / K
    }
    se_cl[s]  <- sqrt(sum(IFc_mrd[, s]^2))
    se_row[s] <- sqrt(sum((fit$IFrow[, , s] %*% g)^2))
  }

  # Multiplier (wild) bootstrap of the linearised estimating function.
  draw_xi <- switch(multiplier,
    rademacher = function(K) sample(c(-1, 1), K, replace = TRUE),
    normal     = function(K) stats::rnorm(K),
    poisson    = function(K) stats::rpois(K, 1) - 1)
  Mstar <- matrix(0, B, nt)
  for (b in seq_len(B)) Mstar[b, ] <- mrd + as.vector(crossprod(IFc_mrd, draw_xi(K)))
  ci_lo <- apply(Mstar, 2, stats::quantile, probs = a, names = FALSE)
  ci_hi <- apply(Mstar, 2, stats::quantile, probs = 1 - a, names = FALSE)

  # Simultaneous band (sup-t): critical value from the max standardised deviation.
  se_pos <- ifelse(se_cl > 0, se_cl, NA)
  Tb <- apply(abs(sweep(Mstar, 2, mrd)) / matrix(se_pos, B, nt, byrow = TRUE), 1,
              max, na.rm = TRUE)
  calpha <- stats::quantile(Tb, level, names = FALSE)

  mk <- function(type, se, lo, hi)
    data.frame(tau = t.hor, target = paste0("mrd_t", t.hor), type = type,
               estimate = mrd, se = se, lower = lo, upper = hi,
               stringsAsFactors = FALSE)

  mrd_df <- rbind(
    mk("analytic_row",     se_row, mrd - z * se_row, mrd + z * se_row),
    mk("analytic_cluster", se_cl,  mrd - z * se_cl,  mrd + z * se_cl),
    mk("multiplier",       apply(Mstar, 2, stats::sd), ci_lo, ci_hi))

  band_df <- mk("band", se_cl, mrd - calpha * se_cl, mrd + calpha * se_cl)
  attr(band_df, "calpha") <- calpha

  list(point = list(t.hor = t.hor, mrd = mrd), mrd = mrd_df, band = band_df)
}
