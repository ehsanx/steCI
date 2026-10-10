#' Analytic standard errors for the constant hazard difference (M1-M3)
#'
#' For the scalar constant-hazard-difference estimand (`const(A.baseline)`),
#' returns the model-based SE, the row-level robust (sandwich) SE, and the
#' cluster-robust SE (clustered on subject `id`) directly from
#' [timereg::aalen()], with Wald confidence intervals.
#'
#' The row-level robust SE is timereg's default and ignores that each subject
#' contributes to several stacked trials; the cluster-robust variant (M3)
#' clusters on the subject. Analytic SEs for the marginal-risk-difference
#' curve require the influence-function/delta method and are provided by
#' [ci_lef()]; the resampling CIs ([ci_bootstrap()], [ci_jackknife()]) also
#' cover the curve.
#'
#' @param dat.seq Output of [ste_prepare()] (must contain `id`).
#' @param weight_col Weight column name (default "ipw.stab.cum").
#' @param level Confidence level (default 0.95).
#' @param ties `"timereg"` (default) fits with [timereg::aalen()], which
#'   breaks tied event times by random jitter; `"simultaneous"` computes the
#'   same estimate and the same three variances (timereg's `var.gamma` and
#'   `robvar.gamma`, including its `max.clust = 1000` cluster binning) with
#'   tied events processed together, via [.semiaalen()]. Identical to machine
#'   precision on data without tied event times. See [ste_fit()].
#'
#' @return A data.frame with columns `method` (model/robust/cluster),
#'   `estimate`, `se`, `lower`, `upper`.
#' @export
ci_analytic <- function(dat.seq, weight_col = "ipw.stab.cum", level = 0.95,
                        ties = c("timereg", "simultaneous")) {
  ties <- match.arg(ties)
  w <- dat.seq[[weight_col]]
  z <- stats::qnorm(1 - (1 - level) / 2)

  if (ties == "simultaneous") {
    sf <- .semiaalen(dat.seq, weight_col = weight_col, variance = TRUE,
                     clusters = dat.seq$id)
    gamma <- sf$gamma[1]
    se <- sqrt(c(sf$var_model[1, 1], sf$var_robust[1, 1], sf$var_cluster[1, 1]))
    return(data.frame(
      method = c("model", "robust", "cluster"), estimate = gamma, se = se,
      lower = gamma - z * se, upper = gamma + z * se, stringsAsFactors = FALSE))
  }

  f0 <- timereg::aalen(
    survival::Surv(time.new, time.stop.new, event) ~
      const(A.baseline) + L.baseline,
    data = dat.seq, n.sim = 0, weights = w, robust = 1)

  gamma <- unname(f0$gamma[1, 1])
  se_model  <- sqrt(f0$var.gamma[1, 1])
  se_robust <- sqrt(f0$robvar.gamma[1, 1])

  # Cluster-robust: refit passing subject id as the cluster variable.
  #
  # max.clust IS DELIBERATELY LEFT AT timereg's DEFAULT. Above a size threshold timereg
  # quantile-bins clusters into a fixed number of groups, so at n = 5000 distinct subjects share
  # bins. It is tempting to pin max.clust = NULL so the code matches the paper's phrase
  # "clustered on the subject" -- but MEASURED on the canonical replicate that changes M3 at
  # n = 5000 from 0.00776678 to 0.00764008, a 1.6% shift, and it is the DEFAULT that reproduces
  # every published value. Pinning would desynchronise the code from the paper.
  #
  # What the binning costs, precisely. The merged meat is
  #   sum_b (sum_{i in b} U_i)(sum_{i in b} U_i)' = sum_i U_i U_i' + sum_b sum_{i != j in b} U_i U_j'.
  # Subjects are independent, so the cross-product term has mean zero and the binned meat is
  # UNBIASED for the exact subject-level meat. It is NOT conservative -- that claim was wrong;
  # on the SE scale Jensen makes it marginally ANTI-conservative. The cost is extra noise, not
  # bias. Note also that timereg's binning test is a strict `<`, so ONLY n = 5000 is binned:
  # n = 300 and n = 1000 are exact subject clusters.
  #
  # MEASURED over 20 replicates at n = 5000 (not one draw -- the single canonical replicate's
  # 1.6% is a one-SD draw and must not be quoted as the effect size): mean relative shift from
  # pinning max.clust = NULL is +0.37%, sd 1.63%, t = 1.01 on 19 df, p = 0.33, 8/20 positive.
  # That is noise about zero, consistent with the unbiasedness above. Mean |shift| is 1.24%,
  # so any SINGLE replicate comparison will look large in either direction -- do not read a
  # systematic effect off one fit. The induced coverage shift is second order and well below
  # the Monte-Carlo SE of the reported entry, so the published 0.855 stands.
  # If a future run wants the unbinned version, pin max.clust = NULL AND re-run -- never pin
  # it alone, or the code will stop reproducing the published table.
  f1 <- timereg::aalen(
    survival::Surv(time.new, time.stop.new, event) ~
      const(A.baseline) + L.baseline,
    data = dat.seq, n.sim = 0, weights = w, robust = 1,
    clusters = dat.seq$id)
  se_cluster <- sqrt(f1$robvar.gamma[1, 1])

  se <- c(se_model, se_robust, se_cluster)
  data.frame(
    method = c("model", "robust", "cluster"),
    estimate = gamma,
    se = se,
    lower = gamma - z * se,
    upper = gamma + z * se,
    stringsAsFactors = FALSE)
}
