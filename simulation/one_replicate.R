# One simulation replicate: simulate a dataset from a scenario, run the point
# estimator and all confidence-interval methods, and return tidy rows with
# schema: estimand, tau, method, estimate, se, lower, upper.
#
# Assumes the steCI package functions are already loaded (via library(steCI)
# or by sourcing steCI/R/*.R). Depends on scenarios.R (DGM_FIXED, T_HOR).

# PER-METHOD FAILURE ISOLATION (fixed 2026-07-17).
#
# THE BUG. Every method block below used to run bare, so ANY of them throwing aborted the whole
# function -- and `run_sim.R`'s single tryCatch then discarded the replicate for EVERY method. Since
# `.aalen_if` (M6) goes singular when a late risk set has no A.baseline variation while
# `timereg::aalen` (M1-M4) survives, the deletion was CONDITIONED ON "M6 did not go singular": a
# non-ignorable mechanism, and it biased the results toward our own recommended method.
#
# MEASURED (`n300_evhigh_cfweak_pvhigh`): 195 of 207 dropped replicates fit M1-M3 fine, i.e.
# the point estimator was constructible in ~98.8% of the cell's replicates.
# Dropped reps had mean gamma -0.0874 vs -0.0362 retained (~1 SD apart); the cell's reported bias was
# +0.0173 against +0.0072 on all reps -- a factor of 2.4.
#
# THE FIX. Each method gets its own tryCatch and its own failure flag. Only `ste_estimate` is a genuine
# prerequisite (M1-M3 need its `ds`, and without a point estimate the replicate has nothing to report);
# `ci_bootstrap`, `ci_lef` and `ci_jackknife` each take `d` directly and are therefore isolated. An M6
# singularity now costs M6's rows ONLY. `failed_methods` records what was lost, so per-method
# failure rates are reported directly rather than inferred.
run_one <- function(cell, rep_seed, t.hor = T_HOR, B = 500,
                    weight = "stab", do_jackknife = FALSE, jack_max_n = 500) {
  set.seed(rep_seed)
  d <- dat_sim(n = cell$n, n.visit = DGM_FIXED$n.visit,
               gamma.0 = cell$gamma.0, gamma.L = cell$gamma.L,
               alpha.0 = cell$alpha.0, alpha.A = DGM_FIXED$alpha.A,
               alpha.L = DGM_FIXED$alpha.L, alpha.U = DGM_FIXED$alpha.U)

  rows <- list()
  add <- function(estimand, tau, method, estimate, se, lower, upper) {
    rows[[length(rows) + 1]] <<- data.frame(
      estimand = estimand, tau = tau, method = method,
      estimate = estimate, se = se, lower = lower, upper = upper,
      stringsAsFactors = FALSE)
  }
  fails <- character(0)
  attempt <- function(name, expr) {
    tryCatch(expr, error = function(e) {
      fails <<- c(fails, name)
      message("  [rep-method] ", name, " failed: ", conditionMessage(e))
      NULL
    })
  }

  # Point estimate + prepared data (reused by the analytic CIs).
  # THE ONLY UNAVOIDABLE DEPENDENCY: no point estimate => the replicate has nothing to contribute.
  est <- attempt("point", ste_estimate(d, n.visit = DGM_FIXED$n.visit, t.hor = t.hor,
                                       weight = weight, estimand = "both"))
  if (is.null(est)) return(NULL)
  ds <- attr(est, "dat.seq")
  add("hd", NA_real_, "point", est$gamma, NA, NA, NA)
  for (j in seq_along(t.hor))
    add("mrd", t.hor[j], "point", est$mrd[j], NA, NA, NA)

  # M1-M3 analytic (hazard difference only).
  # `ds` carries whichever weight column ste_estimate was told to use, so M1-M3 inherit the
  # caller's `weight` argument rather than silently reverting to the untruncated one. The
  # shipped runs are weight = "stab" throughout, so this is a guard on future callers: under
  # weight = "trunc" a bare call would have analysed truncated data with untruncated weights.
  stopifnot(identical(attr(ds, "weight"), weight) || is.null(attr(ds, "weight")))
  an <- attempt("ci_analytic", ci_analytic(ds))
  if (!is.null(an)) for (i in seq_len(nrow(an)))
    add("hd", NA_real_, an$method[i], an$estimate[i], an$se[i], an$lower[i], an$upper[i])

  # M4 nonparametric bootstrap (hd + mrd).
  bo <- attempt("ci_bootstrap", ci_bootstrap(d, n.visit = DGM_FIXED$n.visit, t.hor = t.hor,
                                             weight = weight, estimand = "both", B = B))
  if (!is.null(bo)) {
    for (i in seq_len(nrow(bo$gamma)))
      add("hd", NA_real_, paste0("boot_", bo$gamma$type[i]),
          bo$gamma$estimate[i], bo$gamma$se[i], bo$gamma$lower[i], bo$gamma$upper[i])
    for (i in seq_len(nrow(bo$mrd)))
      add("mrd", bo$mrd$tau[i], paste0("boot_", bo$mrd$type[i]),
          bo$mrd$estimate[i], bo$mrd$se[i], bo$mrd$lower[i], bo$mrd$upper[i])
  }

  # M6 analytic delta-method + multiplier/LEF bootstrap + band (MRD only).
  # This is the block that used to take every other method down with it.
  lf <- attempt("ci_lef", ci_lef(d, n.visit = DGM_FIXED$n.visit, t.hor = t.hor,
                                 weight = weight, B = B))
  if (!is.null(lf)) {
    for (i in seq_len(nrow(lf$mrd)))
      add("mrd", lf$mrd$tau[i], paste0("m6_", lf$mrd$type[i]),
          lf$mrd$estimate[i], lf$mrd$se[i], lf$mrd$lower[i], lf$mrd$upper[i])
    for (i in seq_len(nrow(lf$band)))
      add("mrd", lf$band$tau[i], "m6_band",
          lf$band$estimate[i], lf$band$se[i], lf$band$lower[i], lf$band$upper[i])
  }

  # M6 Approach-2: propagate IPACW weight-estimation uncertainty.
  lf2 <- attempt("ci_lef_a2", ci_lef(d, n.visit = DGM_FIXED$n.visit, t.hor = t.hor,
                                     weight = weight, B = B, approach = "a2"))
  if (!is.null(lf2)) {
    a2 <- lf2$mrd[lf2$mrd$type %in% c("analytic_cluster", "multiplier"), ]
    for (i in seq_len(nrow(a2)))
      add("mrd", a2$tau[i], paste0("m6a2_", a2$type[i]),
          a2$estimate[i], a2$se[i], a2$lower[i], a2$upper[i])
  }

  # M5 jackknife (small n only).
  if (do_jackknife) {
    jk <- attempt("ci_jackknife", ci_jackknife(d, n.visit = DGM_FIXED$n.visit, t.hor = t.hor,
                                               weight = weight, estimand = "both",
                                               max_n = jack_max_n))
    if (!is.null(jk) && !isTRUE(jk$skipped)) {
      add("hd", NA_real_, "jackknife", jk$gamma$estimate, jk$gamma$se,
          jk$gamma$lower, jk$gamma$upper)
      for (i in seq_len(nrow(jk$mrd)))
        add("mrd", jk$mrd$tau[i], "jackknife", jk$mrd$estimate[i], jk$mrd$se[i],
            jk$mrd$lower[i], jk$mrd$upper[i])
    }
  }

  out <- do.call(rbind, rows)
  out$neg_haz_frac <- attr(d, "neg_haz_frac")
  out$boot_n_fail <- if (!is.null(bo)) bo$n_fail else NA_integer_
  # Per-replicate record of what failed: "" when everything ran. This is what makes the per-method
  # failure rates reportable, and what stops a future reader inferring failures from missing rows.
  out$failed_methods <- if (length(fails)) paste(fails, collapse = ";") else ""
  out
}
