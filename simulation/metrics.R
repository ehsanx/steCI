# ADEMP performance measures (Morris, White & Crowther 2019).
# Input: a data.frame of stacked replicate rows (from run_one, with a `rep`
# column and scenario columns) plus a truth lookup.

# truth_lookup: function(estimand, tau, alpha.0) -> true value.
# Compute the per-alpha0 truths (KM survival + marginal hazard difference). This
# is the slow step (aalen on a large randomized sample), so it is cached; see
# aggregate.R. n_pseudo drives the KM survival truth, n_hd (much smaller) the
# marginal-hazard-difference truth.
#
# WARNING (2026-07-19). SUPERSEDED -- do not score published numbers against this path.
# It returns a SINGLE Monte-Carlo draw at a fixed default seed (20222022). Its gamma_true
# is displaced from the exact target by +0.0012 / +0.0023 / +0.0031 at alpha.0 =
# 0.15 / 0.25 / 0.40 -- about one single-draw SD each (0.0015 / 0.0023 / 0.0034), and a
# large fraction of the +0.0073 to +0.0087 E1 bias the study measures.
# The E2/KM truths are displaced too, by up to +0.0016. Being shared across all replicates,
# the displacement is a FIXED offset that does not average out.
#
# Both estimands' truths are available in CLOSED FORM. Use dev/truth_cache_exact.rds for
# anything the paper reports; see dev/truth_closed_form.R for the derivation and
# dev/verify_closed_form.R for the independent checks.
#
# Historical Monte Carlo caches (including data/truth_cache.rds) hold the old
# single draw, exact = FALSE, gamma -0.06064 / -0.05733 / -0.05351.
# aggregate.R now defaults to dev/truth_cache_exact.rds; select a historical
# cache explicitly only when reproducing a diagnostic comparison.
#
# The aalen const() fit itself is not the problem: 2000 independent draws per arm at
# n_hd = 3e4 bound its finite-sample bias below 0.00022 (95%, worst arm). That is a BOUND,
# not a demonstration of unbiasedness, and n_hd was NOT swept -- it says nothing about
# other n_hd. The displacement above is single-draw Monte-Carlo error, not estimator bias.
# This path is retained only so the pipeline is self-contained.
compute_truths <- function(alpha0_levels = c(0.15, 0.25, 0.4), n_pseudo = 5e5,
                           t.hor = 1:5, seed = 20222022, n_hd = 3e4) {
  tv <- lapply(alpha0_levels, function(a0)
    true_values(n = n_pseudo, n.visit = DGM_FIXED$n.visit, t.hor = t.hor,
                alpha.0 = a0, alpha.A = DGM_FIXED$alpha.A,
                alpha.L = DGM_FIXED$alpha.L, alpha.U = DGM_FIXED$alpha.U,
                seed = seed, n_hd = n_hd))
  names(tv) <- as.character(alpha0_levels)
  tv
}

# Build the truth lookup from a computed/cached `tv` list. E1 truth = causal
# MARGINAL constant hazard difference (not the DGM conditional alpha.A: sustained
# treatment also lowers the confounder L).
truth_lookup <- function(tv) {
  function(estimand, tau, alpha.0) {
    t <- tv[[as.character(alpha.0)]]
    if (estimand == "hd") return(t$gamma_true)
    t$mrd[match(tau, t$t.hor)]
  }
}

make_truth <- function(...) truth_lookup(compute_truths(...))

# ---------------------------------------------------------------------------
# Definitional validity of a replicate's outputs.
#
# WHY. A weighted-Aalen fit in a low-overlap cell can emit estimates OUTSIDE the
# parameter space and near-infinite SEs. Measured on the 2026-07 run: 36 mrd point
# estimates violate |MRD| <= 1 (to -2.587), and bootstrap SEs reach 1.48e6 against a
# typical 0.07. Left unguarded these poison the SCALE statistics: emp_sd is sd() over
# them, so every se_ratio inherits the contamination (one replicate inflated a cell's
# emp_sd by 16%).
#
# MRD = Pr(T^{a=1}<=tau) - Pr(T^{a=0}<=tau) is a difference of probabilities, so it
# lives in [-1,1] BY DEFINITION -- an estimate outside it is an estimator failure, not
# data. Its SE cannot legitimately exceed `se_max_mrd`: a 95% interval of half-width
# 1.96*0.5 already spans the entire parameter space.
#
# The hazard difference (estimand "hd") has no such definitional bound, so only
# finiteness is required there.
#
# These are GUARDS ON SCALE STATISTICS, NOT A TRIM. Coverage is deliberately left
# unguarded below: dropping replicates that sit far from the centre removes exactly the
# ones that fail to cover, which inflates coverage even on perfectly Gaussian data.
# Invalid replicates are COUNTED and reported, never silently discarded.
# ---------------------------------------------------------------------------
.valid_estimate <- function(estimand, estimate) {
  is.finite(estimate) & (estimand != "mrd" | abs(estimate) <= 1)
}
.valid_se <- function(estimand, se, se_max_mrd = 0.5) {
  is.finite(se) & se > 0 & (estimand != "mrd" | se <= se_max_mrd)
}

# Coverage-relevant scale: for a symmetric interval of half-width h recentred on
# mean_est, coverage_be = 0.95 exactly when h = quantile(|est - mean_est|, 0.95).
# Dividing by qnorm(0.975) puts it on the SD scale (it equals sd() for a Gaussian),
# so mean_se/fold_sd is directly comparable to se_ratio but is not distorted by a
# heavy tail. Chosen as the scale that reproduces measured coverage
# (mean_se/emp_sd and mean_se/mad_sd do not).
.fold_sd <- function(est, centre) {
  if (length(est) < 2L) return(NA_real_)
  unname(stats::quantile(abs(est - centre), 0.95, na.rm = TRUE) / stats::qnorm(0.975))
}

# Compute metrics per (scenario, estimand, tau, method).
summarize_metrics <- function(reps, truth_fun, se_max_mrd = 0.5) {
  reps$truth <- mapply(truth_fun, reps$estimand, reps$tau, reps$alpha.0)

  # Empirical bias / SD from the point estimates.
  pt <- reps[reps$method == "point", ]
  key <- function(df) paste(df$scenario, df$estimand, df$tau, sep = "|")
  pt_split <- split(pt, key(pt))
  emp <- do.call(rbind, lapply(pt_split, function(g) {
    fin <- g$estimate[is.finite(g$estimate)]           # legacy basis (contaminated)
    ok  <- .valid_estimate(g$estimand, g$estimate)
    est <- g$estimate[ok]                              # guarded basis
    mean_est <- mean(est)
    data.frame(scenario = g$scenario[1], estimand = g$estimand[1], tau = g$tau[1],
               truth = g$truth[1],
               n_sim = length(fin),
               n_valid = length(est),
               n_invalid_est = sum(is.finite(g$estimate)) - length(est),
               mean_est = mean_est, bias = mean_est - g$truth[1],
               emp_sd = stats::sd(est),
               fold_sd = .fold_sd(est, mean_est),
               mse = mean((est - g$truth[1])^2),
               # unguarded values, retained so the contamination stays auditable
               mean_est_raw = mean(fin), emp_sd_raw = stats::sd(fin),
               stringsAsFactors = FALSE)
  }))
  emp_key <- paste(emp$scenario, emp$estimand, emp$tau, sep = "|")

  # CI-method metrics.
  ci <- reps[reps$method != "point", ]
  ci_split <- split(ci, paste(ci$scenario, ci$estimand, ci$tau, ci$method, sep = "|"))
  out <- do.call(rbind, lapply(ci_split, function(g) {
    truth <- g$truth[1]
    ek <- paste(g$scenario[1], g$estimand[1], g$tau[1], sep = "|")
    i <- match(ek, emp_key)
    esd <- emp$emp_sd[i]; fsd <- emp$fold_sd[i]
    esd_raw <- emp$emp_sd_raw[i]; mean_est <- emp$mean_est[i]
    # Coverage is NOT guarded -- see the note above .valid_estimate().
    covered <- g$lower <= truth & g$upper >= truth
    covered <- covered[is.finite(covered)]
    covered_be <- g$lower <= mean_est & g$upper >= mean_est
    covered_be <- covered_be[is.finite(covered_be)]
    ns <- length(covered)
    cov <- mean(covered)
    se_ok <- .valid_se(g$estimand, g$se, se_max_mrd)
    se_v <- g$se[se_ok]
    mean_se <- mean(se_v)
    mean_se_raw <- mean(g$se[is.finite(g$se)])
    data.frame(
      scenario = g$scenario[1], estimand = g$estimand[1], tau = g$tau[1],
      method = g$method[1], truth = truth, n_sim = ns,
      coverage = cov,
      coverage_be = mean(covered_be),
      mc_se_cov = sqrt(cov * (1 - cov) / ns),
      mean_se = mean_se, median_se = stats::median(se_v),
      emp_sd = esd, fold_sd = fsd,
      se_ratio = mean_se / esd,
      se_ratio_fold = mean_se / fsd,
      n_invalid_se = sum(is.finite(g$se)) - length(se_v),
      n_invalid_est = emp$n_invalid_est[i],
      # unguarded values, retained so the contamination stays auditable
      mean_se_raw = mean_se_raw, emp_sd_raw = esd_raw,
      se_ratio_raw = mean_se_raw / esd_raw,
      mean_width = mean((g$upper - g$lower)[is.finite(g$upper - g$lower)]),
      stringsAsFactors = FALSE)
  }))

  list(empirical = emp, ci = out, joint = .joint_coverage(reps))
}

# Joint (simultaneous) coverage across horizons: the fraction of replicates whose
# interval covers the truth at EVERY tau at once. The paper's band claim is a
# simultaneous claim, but only pointwise coverage was ever computed -- we measured
# the true value at 0.8807 against a published pointwise 0.960. A curve counts only if
# all horizons are present, so partial curves cannot inflate it.
.joint_coverage <- function(reps) {
  r <- reps[reps$method != "point" & !is.na(reps$tau) &
            is.finite(reps$lower) & is.finite(reps$upper), ]
  if (!nrow(r)) return(NULL)
  n_tau <- length(unique(r$tau))
  if (n_tau < 2L) return(NULL)
  r$cov <- r$lower <= r$truth & r$upper >= r$truth
  per_rep <- stats::aggregate(list(k = r$cov, all_cov = r$cov),
                              by = list(scenario = r$scenario, estimand = r$estimand,
                                        method = r$method, rep = r$rep),
                              FUN = function(x) c(length(x), all(x)))
  k <- per_rep$k[, 1]; allc <- per_rep$all_cov[, 2]
  keep <- k == n_tau
  d <- data.frame(scenario = per_rep$scenario[keep], estimand = per_rep$estimand[keep],
                  method = per_rep$method[keep], all_cov = as.logical(allc[keep]),
                  stringsAsFactors = FALSE)
  o <- stats::aggregate(all_cov ~ scenario + estimand + method, data = d,
                        FUN = function(x) c(joint_coverage = mean(x), n_curves = length(x)))
  data.frame(o[, c("scenario", "estimand", "method")],
             joint_coverage = o$all_cov[, 1], n_curves = o$all_cov[, 2],
             mc_se = sqrt(o$all_cov[, 1] * (1 - o$all_cov[, 1]) / o$all_cov[, 2]),
             stringsAsFactors = FALSE)
}
