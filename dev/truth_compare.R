# =====================================================================
# truth_compare.R -- compare a CLOSED-FORM E1 target g_star against the
# ARC array's Monte-Carlo estimates, and in doing so MEASURE the
# timereg::aalen const() estimator's own finite-sample bias in n_hd.
#
#   Rscript dev/truth_collect.R      # writes truth_draws_all.rds
#   Rscript dev/truth_compare.R      # this file
#
# WHAT IS BEING COMPARED, PRECISELY.
# The array does not estimate "the population gamma". It estimates
#   E[ gamma_hat(n_hd) ],
# the expectation of the timereg estimator AT a finite n_hd. The closed form
# g_star is the n -> Inf limit of that same estimator. They differ by the
# estimator's finite-sample bias, which is O(1/n) because gamma_hat is a ratio
# of two sample averages. So the correct comparison is NOT g_star against any
# single n_hd arm; it is g_star against the INTERCEPT of a fit of
#   g_bar(n_hd)  ~  g_star + b / n_hd
# across the sweep. The slope b is not a nuisance -- it is a reportable result.
#
# g_star MUST be the limit of the SAME estimating equation the array runs.
# See the header of dev/truth_closed_form.R for that derivation; the short
# version is that for aalen(Surv(Tstop,event) ~ const(A)) on a 50/50 sample,
#   g_star = int_0^tau w(t) [lambda_1(t) - lambda_0(t)] dt / int_0^tau w(t) dt,
#   w(t) = S_1(t) S_0(t) / (S_1(t) + S_0(t)).
# A g_star computed as anything else -- e.g. an unweighted time-average of the
# hazard difference, or alpha.A -- will disagree with the array, and the
# disagreement will be MISREAD as aalen bias. Guard that with the S_a(t)
# prerequisite check (check_survival_first()) before trusting any gamma result.
# =====================================================================

# ---------------------------------------------------------------------
# GLS fit of y = g0 + b*x with KNOWN variances s^2 (not lm(): lm rescales the
# SEs by the residual mean square, which is wrong when s is known a priori --
# and here s is known, it is the Monte-Carlo SE from K draws).
# ---------------------------------------------------------------------
.gls_line <- function(x, y, s) {
  X <- cbind(1, x); W <- 1 / s^2
  XtWX <- crossprod(X, X * W)
  V    <- solve(XtWX)
  beta <- drop(V %*% crossprod(X, y * W))
  res  <- y - drop(X %*% beta)
  list(intercept = beta[1], slope = beta[2],
       se_intercept = sqrt(V[1, 1]), se_slope = sqrt(V[2, 2]),
       chisq = sum(res^2 / s^2), df = length(y) - 2L, resid = res)
}

#' Compare a closed-form E1 target against the ARC truth array.
#'
#' @param d        data.frame with columns alpha0, n_hd, gamma -- i.e. exactly
#'                 what dev/truth_collect.R saves to truth_draws_all.rds.
#' @param g_star   named numeric, closed-form target per alpha.0. Names must be
#'                 the alpha.0 values as character, e.g. c("0.15"=..., "0.25"=...).
#' @param quad_err absolute quadrature error bound on g_star (scalar or named).
#'                 Folded into the comparison SE so a sloppy quadrature cannot
#'                 masquerade as agreement. Default 1e-9 (deterministic nested
#'                 quadrature should beat this easily; if it does not, say so).
#' @param e1_bias  the E1 |bias| the paper reports, used only to express
#'                 discrepancies on a scale a reader cares about.
#' @param prac_frac practical-equivalence bound as a fraction of e1_bias.
#'                 With K = 2000 the MCSE is ~1e-4, so a discrepancy can be
#'                 statistically significant and scientifically irrelevant.
#'                 Report BOTH verdicts; do not let either alone decide.
truth_compare <- function(d,
                          g_star,
                          quad_err   = 1e-9,
                          e1_bias    = 0.0057,
                          prac_frac  = 0.10,
                          z_crit     = 3) {

  stopifnot(all(c("alpha0", "n_hd", "gamma") %in% names(d)))
  d  <- d[is.finite(d$gamma), ]
  a0s <- sort(unique(d$alpha0))
  hs  <- sort(unique(d$n_hd))
  prac <- prac_frac * e1_bias

  gs <- function(a) {
    v <- g_star[[as.character(a)]]
    if (is.null(v) || !is.finite(v)) stop("no closed-form g_star for alpha.0 = ", a)
    v
  }
  qe <- function(a) if (length(quad_err) > 1) quad_err[[as.character(a)]] else quad_err

  # ---- per-cell summary --------------------------------------------------
  cells <- do.call(rbind, lapply(a0s, function(a) {
    do.call(rbind, lapply(hs, function(h) {
      g <- d$gamma[d$alpha0 == a & d$n_hd == h]
      if (!length(g)) return(NULL)
      m <- mean(g); s <- sd(g) / sqrt(length(g))
      data.frame(alpha0 = a, n_hd = h, K = length(g), mean = m,
                 sd1 = sd(g), mcse = s, g_star = gs(a),
                 bias = m - gs(a),
                 z = (m - gs(a)) / sqrt(s^2 + qe(a)^2),
                 n_times_bias = h * (m - gs(a)),
                 pct_e1 = 100 * (m - gs(a)) / e1_bias)
    }))
  }))

  cat("\n=========== (1) closed form vs array, per (alpha.0, n_hd) ===========\n")
  cat("bias = E[gamma_hat(n_hd)] - g_star. If the O(1/n) model holds, n*bias is\n")
  cat("roughly CONSTANT down each alpha.0 block -- that column is the real check.\n\n")
  cat(sprintf("%-7s %8s %5s %11s %9s %11s %11s %8s %11s %8s\n",
              "alpha.0","n_hd","K","g_bar","MCSE","g_star","bias","bias/SE","n*bias","%E1bias"))
  for (i in seq_len(nrow(cells))) with(cells[i, ],
    cat(sprintf("%-7s %8d %5d %11.6f %9.6f %11.6f %+11.6f %8.1f %+11.2f %+8.1f\n",
                format(alpha0), n_hd, K, mean, mcse, g_star, bias, z, n_times_bias, pct_e1)))

  # ---- extrapolation to n = Inf -----------------------------------------
  cat("\n=========== (2) extrapolated limit vs closed form ===========\n")
  if (length(hs) < 2) {
    cat("  only one n_hd present (the 'quick' profile). The O(1/n) bias CANNOT be\n")
    cat("  separated from a wrong closed form. Treat the row above as an UPPER BOUND\n")
    cat("  on |g_star - limit|, not as a test. Run the 'full' profile.\n")
    return(invisible(list(cells = cells, fits = NULL)))
  }
  cat("GLS fit of g_bar(n_hd) = g0 + b/n_hd with known MC variances.\n")
  cat("Verdict compares the INTERCEPT g0 to g_star -- that is the like-for-like test.\n\n")
  cat(sprintf("%-7s %12s %10s %11s %8s %12s %10s %s\n",
              "alpha.0","g0 (extrap)","SE(g0)","g_star","z","b (slope)","LoF chi2","verdict"))
  fits <- list(); zs <- numeric(0)
  for (a in a0s) {
    cc <- cells[cells$alpha0 == a, ]
    if (nrow(cc) < 2) next
    f  <- .gls_line(1 / cc$n_hd, cc$mean, cc$mcse)
    se <- sqrt(f$se_intercept^2 + qe(a)^2)
    z  <- (f$intercept - gs(a)) / se
    zs <- c(zs, z)
    lof <- if (f$df > 0) sprintf("%.2f(%d)", f$chisq, f$df) else "n/a(0)"
    verdict <- if (abs(z) < z_crit) "AGREE"
               else if (abs(f$intercept - gs(a)) < prac) "stat. sig. but < practical bound"
               else "DISAGREE"
    cat(sprintf("%-7s %12.6f %10.6f %11.6f %+8.2f %+12.1f %10s  %s\n",
                format(a), f$intercept, se, gs(a), z, f$slope, lof, verdict))
    fits[[as.character(a)]] <- f
  }

  # Joint test: three alpha.0 levels are independent runs, so z^2 sums to chi2_3.
  if (length(zs)) {
    X2 <- sum(zs^2); p <- pchisq(X2, length(zs), lower.tail = FALSE)
    cat(sprintf("\n  JOINT over %d alpha.0 levels: chi2 = %.2f on %d df, p = %.3g  -> %s\n",
                length(zs), X2, length(zs), p,
                if (p > 0.001) "consistent with the closed form" else "NOT consistent"))
  }

  # ---- assumption-light fallback ----------------------------------------
  cat("\n=========== (3) assumption-light cross-check ===========\n")
  cat("If the lack-of-fit chi2 above is large, the 1/n model is wrong and the\n")
  cat("extrapolation is not trustworthy. Then fall back to this, which assumes\n")
  cat("only that |bias| is DECREASING in n_hd -- so the largest-n_hd gap bounds it.\n\n")
  big <- max(hs)
  for (a in a0s) {
    cc <- cells[cells$alpha0 == a & cells$n_hd == big, ]
    if (!nrow(cc)) next
    cat(sprintf("  alpha.0=%-5s |g_bar(%d) - g_star| = %.6f  (%.1f MCSE, %.1f%% of the E1 bias)  %s\n",
                format(a), big, abs(cc$bias), abs(cc$z), abs(cc$pct_e1),
                if (abs(cc$bias) < prac) "within practical bound" else "exceeds practical bound"))
  }

  cat(sprintf("\n  practical-equivalence bound = %.1f%% of E1 |bias| %.4f = %.6f\n",
              100 * prac_frac, e1_bias, prac))
  cat("\n=========== (4) what to ship ===========\n")
  cat("If (2) AGREES: ship g_star in the truth cache -- it is exact, and the array's\n")
  cat("  residual MCSE then contributes ZERO to the reported E1 bias.\n")
  cat("If (2) DISAGREES: do NOT ship either number yet. A disagreement is far more\n")
  cat("  likely to be a wrong WEIGHT FUNCTION in the derivation than a wrong array --\n")
  cat("  run check_survival_first() to localise it before touching the cache.\n")

  invisible(list(cells = cells, fits = fits, prac = prac))
}


#' PREREQUISITE. Localise a disagreement before interpreting it.
#'
#' The closed form has two independent halves that can each be wrong:
#'   (A) the DGM half  -- S_a(t) and f_a(t), the marginal survival and density
#'       under a constant arm, obtained from the Gaussian L-recursion with the
#'       hazard replaced by its positive part (see note below);
#'   (B) the ESTIMATOR half -- the weight w(t) = S_1 S_0/(S_1+S_0) that turns
#'       those into the probability limit of timereg's estimating equation.
#'
#' S_a(t) depends ONLY on (A). The paper already has high-precision KM estimates
#' of S_a at t = 1..5 (n = 1e6 per arm -> SE ~ 5e-4). So:
#'   S_a matches, gamma does not  ->  the bug is in (B), the weight function.
#'   S_a does not match           ->  the bug is in (A); gamma cannot be right either.
#' Run this FIRST. Without it, a gamma disagreement is uninterpretable.
#'
#' @param s_star  closed-form S_a(t) at t.hor, list(surv1=, surv0=) per alpha.0
#' @param s_mc    the Monte-Carlo KM values, same shape
#' @param se_mc   their SEs (KM SE at n per arm; ~5e-4 at n = 1e6)
check_survival_first <- function(s_star, s_mc, se_mc, t.hor = 1:5, z_crit = 3) {
  cat("\n=========== PREREQUISITE: does the DGM half reproduce S_a(t)? ===========\n")
  worst <- 0
  for (a in names(s_star)) for (arm in c("surv1", "surv0")) {
    ds <- s_star[[a]][[arm]] - s_mc[[a]][[arm]]
    z  <- ds / se_mc[[a]][[arm]]
    worst <- max(worst, max(abs(z)))
    cat(sprintf("  alpha.0=%-5s %s  max|diff| = %.6f  max|z| = %.2f\n",
                a, arm, max(abs(ds)), max(abs(z))))
  }
  cat(sprintf("\n  worst |z| over all (alpha.0, arm, t) = %.2f -> %s\n", worst,
              if (worst < z_crit) "DGM half VERIFIED; any gamma gap is the WEIGHT function"
              else "DGM half FAILS; fix it before looking at gamma at all"))
  invisible(worst)
}


# ---------------------------------------------------------------------
if (sys.nframe() == 0L) {
  getc <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) x else d }
  ROOT <- getc("REPO_ROOT", ".")
  OUT  <- getc("TRUTH_OUT", file.path(ROOT, "simulation", "output", "truth"))
  f    <- file.path(OUT, "truth_draws_all.rds")
  if (!file.exists(f)) stop("run dev/truth_collect.R first; no ", f)
  d <- readRDS(f)

  gf <- getc("TRUTH_GSTAR", file.path(ROOT, "dev", "gamma_closed_form.rds"))
  if (!file.exists(gf))
    stop("no closed-form values at ", gf, "\n",
         "  Expected a named numeric like c(`0.15`=..., `0.25`=..., `0.4`=...),\n",
         "  computed as the probability limit of aalen(Surv ~ const(A)) -- NOT as a\n",
         "  plain time-average of the hazard difference. See this file's header.")
  g_star <- readRDS(gf)

  truth_compare(d, g_star)
}
