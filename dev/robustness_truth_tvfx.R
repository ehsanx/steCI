# EXACT-TO-TOLERANCE truth for robustness DGM (ii): the TIME-VARYING TREATMENT EFFECT
# generator (dev/robustness_tvfx.R, Web Table S7).  Replaces the single-frozen-draw
# pseudo-RCT truth `true_gamma()` (n=3e4/arm, seed 2024) with deterministic quadrature.
#
# This is the third and last of the three robustness truths, and the only one where the
# `h > 0` guard is a FIRST-ORDER feature of the data-generating process rather than a
# correction term.  It is therefore handled EXACTLY, not bounded.
#
# ---------------------------------------------------------------------------------------
# (A) WHAT THE GUARD DOES -- read off dev/robustness_tvfx.R:31-36, not assumed.
#
#       for (k in 1:NV) { u <- runif(n); h <- tv_haz(k, A, L, U, alpha.0); nt <- -log(u)/h
#         Tt <- ifelse(is.na(Tt) & nt < 1 & h > 0, k - 1 + nt, Tt) }
#
#     When h_k > 0, nt = -log(u)/h_k and the event lands in [k-1,k) iff nt < 1 iff
#     u > exp(-h_k): the interval contributes survival factor exp(-h_k).
#     When h_k <= 0, nt = -log(u)/h_k <= 0, so `nt < 1` is satisfied but `h > 0` is NOT,
#     and the row is skipped: NO event can occur in that interval, and the subject stays
#     at risk.  The interval contributes survival factor 1 = exp(-0).
#     Both branches are exp(-max(h_k, 0)).  So the process actually simulated is the one
#     with hazard  h_k^+ = max(h_k, 0)  held constant on [k-1, k), and
#
#         S_a(t) = E[ exp( - sum_k w_k * max(h_k, 0) ) ],   w = interval occupancy up to t.
#
#     *** This is an IDENTITY for the generator, not an approximation. ***  There is no
#     guard error term to bound here, and no closed form to correct: the truncation is
#     inside the expectation and is integrated exactly.  (Contrast dev/robustness_truth_exact.R,
#     where the closed form integrates the UNCLAMPED hazard and needs a measured +4.1e-5 CRN
#     correction at alpha0=0.15.  That framing is not merely inconvenient here, it is wrong:
#     max_k P(h_k <= 0) reaches 0.53 in the treated arm at alpha0=0.15 -- see CHECK (c).)
#
#     Administrative censoring: `Ts = ifelse(is.na(Tt), NV, Tt)`, D = 0 there, in BOTH arms.
#     So tau = NV = 5 and K = 5 intervals.
#
# (B) NO CLOSED FORM.  The hazard
#         h_k = alpha0 + (AA + BETA*(k-1)) * a + AL * L_k + AU * U
#     IS additive and linear in a jointly Gaussian vector, so E[exp(-sum w_k h_k)] would be
#     an elementary lognormal MGF exp(-mu + v/2).  But  x |-> max(x,0)  is not linear, so
#     E[exp(-sum w_k h_k^+)] is a genuine K-dimensional Gaussian expectation of a KINKED
#     integrand.  The time-varying coefficient (AA + BETA*(k-1)) is DETERMINISTIC given the
#     arm: it shifts the MEAN of h_k per interval and leaves the covariance untouched.
#
# (C) THE LAW OF (h_1..h_K) UNDER A SUSTAINED ARM a.
#     Under arm a the confounder recursion is
#         L_1 = U + e_1,   L_k = 0.8 L_{k-1} - a + 0.1(k-1) + U + e_k,   e_k ~ N(0,1) iid.
#     Split off the shared frailty, L_k = M_k + c_k U, with
#         c_1 = 1,   c_k = 0.8 c_{k-1} + 1
#         M_1 ~ N(0,1),   M_k = 0.8 M_{k-1} - a + 0.1(k-1) + e_k     <-- U-FREE, and MARKOV
#     so E M_k = b_k (b_1 = 0, b_k = 0.8 b_{k-1} - a + 0.1(k-1)) and Var M_k = vR_k
#     (vR_1 = 1, vR_k = 0.64 vR_{k-1} + 1), Cov(M_j,M_k) = 0.8^|k-j| vR_min(j,k).  Then
#         h_k = A_k + AL * M_k + g_k * U,    A_k = alpha0 + (AA + BETA(k-1)) a,
#                                            g_k = AL * c_k + AU,
#         E h_k    = A_k + AL * b_k
#         Cov(h_j,h_k) = AL^2 * 0.8^|k-j| vR_min(j,k)  +  g_j g_k * SIGU^2.
#
# (D) HOW IT IS EVALUATED, and why not tensor Gauss-Hermite.
#     dev/robustness_truth_cox.R uses tensor GH because its integrand is ENTIRE in every
#     coordinate and nearly affine over the support, so GH converges super-geometrically.
#     Here max(h_k,0) has a KINK in every coordinate, sitting right in the bulk of the mass
#     (P(h_5 <= 0) = 0.53 in the treated arm at alpha0=0.15).  Gauss-Hermite is a polynomial
#     rule; against a kink its convergence degrades from super-geometric to algebraic, and a
#     refinement sequence in m would be slow and non-monotone.  So tensor GH is used ONLY as
#     an independent cross-check (CHECK (d)), not as the production route.
#
#     PRODUCTION ROUTE: the M-chain is Gauss-Markov and -- by construction in (C) -- its
#     transition kernel does NOT depend on U.  U enters only by shifting each interval's
#     kink location, g_k * U.  So:
#        * outer 1-D Gauss-Hermite over U (n_U nodes),
#        * inner FORWARD FILTER over a fine grid in M, propagating the sub-density
#              p_1(l) = phi(l)
#              p_{k+1}(l') = int p_k(l) exp(-h_k(l)^+) N(l'; 0.8 l + (-a + 0.1 k), 1) dl
#          and reading off, for t in interval K with s = t - (K-1),
#              S_a(t) = int p_K(l) exp(-s h_K(l)^+) dl
#              f_a(t) = -dS_a/dt = int p_K(l) h_K(l)^+ exp(-s h_K(l)^+) dl.
#     f is obtained by differentiating UNDER the integral (one extra factor h_K^+ in the same
#     quadrature), exactly as the Cox derivation does -- never by differencing S.  Note S is
#     ANALYTIC in s within an interval even though the integrand is kinked in l: the kink is in
#     the l-variable, which is integrated out, and enters s only through the fixed value h_K^+.
#     Cost: the K-1 transition matrices are built once per arm, and every t inside interval K
#     reuses the same p_K -- so 120 t-nodes cost 4 matrix products, not 120 recursions.
#
#     THE KINK AND THE l-QUADRATURE.  The l-integrand IS kinked (at l = -(A_k + g_k u)/AL),
#     so a high-order rule buys nothing across it; a fine composite trapezoid is used and
#     convergence is DEMONSTRATED by halving the mesh (CHECK (e)).  Honest statement: the
#     kink costs roughly one order of convergence relative to a smooth integrand, which is
#     why nL is in the thousands rather than the dozens.  It is still deterministic, positive
#     weighted, and auditable.
#
# (E) F2 -- gamma is a functional of S_1 and S_0 ALONE.  VERIFIED, verbatim from
#     dev/robustness_truth_cox.R (C).  gamma_true is DEFINED as the population solution of the
#     Lin-Ying constant-additive-effect estimating equation that
#     `aalen(Surv(Ts,ev) ~ const(A), n.sim=0)` solves on the pseudo-RCT.  Profiling out the
#     nonparametric baseline gives
#         0 = int sum_i (A_i - Abar(t)) [dN_i(t) - Y_i(t) A_i gamma dt],
#         Abar(t) = sum_i Y_i(t) A_i / sum_i Y_i(t),
#     and with equal allocation, A randomised / time-fixed / independent of (U,L), and
#     administrative censoring at tau=5 in BOTH arms,
#         Abar(t) -> S1/(S1+S0),   E[dN | A=a] = f_a dt,   E[Y | A=a] = S_a,
#         gamma = int_0^tau [S0 f1 - S1 f0]/(S1+S0) dt  /  int_0^tau S1 S0/(S1+S0) dt.
#     Every step is a property of the ESTIMATOR and of the pseudo-RCT DESIGN.  Nothing refers
#     to how T was generated -- the only covariate is A, so the estimating equation never needs
#     E[covariate | at risk].  In particular the GUARD, which changes only the LAW of T, cannot
#     touch this algebra: it changes S_1 and S_0 and nothing else, and S_1/S_0 are exactly what
#     the functional consumes.  A time-varying true effect is what makes gamma a working-model
#     SUMMARY rather than a structural parameter; it does not alter the functional.
#     F2 HOLDS.  It is not taken on faith: CHECK (f) runs 600 independent pseudo-RCTs through
#     the study's OWN code path and finds |z| < 2 at all three alpha0.
#
# ---------------------------------------------------------------------------------------
# RESULTS (K = 5, tau = 5).  These are what dev/truth_tvfx_exact.rds holds.
#
#     alpha0     gamma_exact      frozen seed-2024 draw    displacement
#     0.15      -0.0954857381         -0.0942775719          +0.001208
#     0.25      -0.0917840129         -0.0898514345          +0.001933
#     0.40      -0.0832047917         -0.0837105519          -0.000506
#
#   The frozen draw's displacement is 0.88 / 0.87 / 0.17 single-draw MC SD -- unremarkable AS A
#   DRAW, which is exactly why it was easy to miss.  It is nevertheless a FIXED offset shared by
#   all 500 replicates of every cell, and it is 21 / 21 / 4 times the MCSE of the 600-seed mean
#   below.  At alpha0=0.25 it is 2.1% of |gamma|.  That is the defect being removed.
#
# ACCURACY.
#   * The guard is handled EXACTLY (section (A)), so there is NO guard error term and NO
#     Monte-Carlo ingredient anywhere in the truth.  Every digit is deterministic quadrature.
#   * Deterministic quadrature, converged and independently confirmed.  The observed convergence
#     floor -- NOT a proven bound -- is:
#       - l-mesh: gamma over nL = 800 / 1600 / 3200 spreads 5e-9 at alpha0=0.15 (the kink limits
#         it) and <= 3e-12 at 0.25 and 0.40.
#       - U-rule: S_1(5) over nU = 16 / 24 / 32 / 48 spreads 1.3e-8 at 0.15, <= 1e-11 elsewhere.
#       - l-range: widening [-16,14] -> [-24,22] moves S_1(5) by <= 5e-9 at 0.15, <= 2e-13 else.
#       - t-rule: gamma over nGL = 8 / 16 / 24 / 32 is identical to 13 decimals at every alpha0.
#     So quote gamma to ~1e-8 at alpha0=0.15 and ~1e-11 at 0.25/0.40.  That is 5 orders of
#     magnitude below the 1.2e-3 displacement it is being used to detect.
#   * INDEPENDENT confirmations, each able to fail:
#       - CHECK (b): where the guard is inactive (P(h_k<=0) <= 5e-27), the filter reproduces the
#         elementary lognormal MGF exp(-E + Var/2) to <= 3.3e-13 absolute at t = 1, 2.5, 5 in
#         both arms.  This exercises the whole 5-step recursion, kernel, l-grid and U-rule while
#         switching the guard OFF, separating "the Gaussian algebra is right" from "the guard
#         handling is right".  At alpha0=0.25 arm 1, where P(h_5<=0)=3.2e-5, the two DIVERGE by
#         -2.4e-7 with the correct sign (max(h,0) >= h => S_guarded <= S_unguarded) -- the guard
#         being seen directly.
#       - CHECK (d): a disjoint reimplementation (tensor Gauss-Hermite on the full 5-D Cov +
#         Cholesky, no Markov split, no l-grid) agrees to <= 1e-9 in the five configurations
#         where the guard is inactive.  In the sixth -- treated arm, alpha0=0.15 -- it agrees
#         only to ~2e-6 and is NON-MONOTONE in m (m=8/12/18 straddle the filter value).  Reported
#         honestly: that is the predicted algebraic-not-spectral convergence of a polynomial rule
#         against a kink, and it is the reason tensor GH is the cross-check and not the route.
#       - f = -dS/dt (differentiated under the integral) vs a central difference of S:
#         <= 2e-12 at t=3.4, all arms and alpha0.
#       - ly_two_sample_tv vs timereg::aalen(~ const(A)) on this DGM at n=4000: measured at
#         4e-16 / 3e-17 / 4e-16.  MEASURED, not claimed as an algebraic identity -- an
#         independent reimplementation of the 2confounder study showed a systematically signed
#         ~1e-6 offset at n=2e4/arm, so this is checked per-DGM and reported, never assumed.
#
# CHECK (a)/(f), the primary check -- 600 independent pseudo-RCTs per alpha0 at n=3e4/arm,
#   generated by sys.source()ing dev/robustness_tvfx.R's OWN gen/tv_haz/event_times (not a
#   retyping) and scored by ly_two_sample_tv.  The driver's `--mc` branch below instead uses the
#   transcribed pseudo_rct_tv(); the two are confirmed to be the SAME generator because their
#   seed-2024 frozen draws agree to every printed digit at all three alpha0
#   (-0.09427757 / -0.08985143 / -0.08371055).  Numbers below are from the sys.source()d run:
#
#     alpha0    MC mean (600 seeds)     MCSE      1-draw sd     exact          z
#     0.15         -0.09541774        5.63e-05     0.001379   -0.09548574    +1.21
#     0.25         -0.09175383        9.02e-05     0.002210   -0.09178401    +0.33
#     0.40         -0.08313668        1.23e-04     0.003010   -0.08320479    +0.55
#
#   All |z| < 1.3 against a truth carrying no MC error of its own.  With 53% of the treated arm's
#   k=5 hazard truncated at alpha0=0.15, this is the check that would have caught a wrong guard
#   characterisation, and it does not fire.  The three alpha0 share seed blocks, so their z are
#   correlated; none is near significance, so no second block was run.
# ---------------------------------------------------------------------------------------
#
#   Rscript dev/robustness_truth_tvfx.R          # derivation + checks (b)-(e), writes the .rds
#   Rscript dev/robustness_truth_tvfx.R --mc     # ... plus CHECK (a)/(f), multi-seed MC (slow)
# ---------------------------------------------------------------------------------------

NV_TV   <- 5L
AA_TV   <- -0.04; BETA_TV <- -0.02; AL_TV <- 0.015; AU_TV <- 0.015
SIGU_TV <- 0.1

# ---- quadrature rules (Golub-Welsch; no package dependency) ---------------------------
gh_rule_tv <- function(m) {          # E[g(X)], X ~ N(0,1); weights sum to 1
  b <- sqrt(seq_len(m - 1) / 2)
  J <- matrix(0, m, m); J[cbind(1:(m-1), 2:m)] <- b; J[cbind(2:m, 1:(m-1))] <- b
  ev <- eigen(J, symmetric = TRUE); o <- order(ev$values)
  list(x = sqrt(2) * ev$values[o], w = ev$vectors[1, o]^2)
}
gl_rule_tv <- function(m) {          # Gauss-Legendre on [-1,1]
  i <- seq_len(m - 1); b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, m, m); J[cbind(1:(m-1), 2:m)] <- b; J[cbind(2:m, 1:(m-1))] <- b
  ev <- eigen(J, symmetric = TRUE); o <- order(ev$values)
  list(x = ev$values[o], w = 2 * ev$vectors[1, o]^2)
}

# ---- section (C): moments under a sustained arm ---------------------------------------
.moments_tv <- function(a) {
  cU <- b <- vR <- numeric(NV_TV); cU[1] <- 1; b[1] <- 0; vR[1] <- 1
  for (k in 2:NV_TV) {
    cU[k] <- 0.8 * cU[k - 1] + 1
    b[k]  <- 0.8 * b[k - 1] - a + 0.1 * (k - 1)
    vR[k] <- 0.64 * vR[k - 1] + 1
  }
  ix <- seq_len(NV_TV)
  list(cU = cU, b = b, vR = vR,
       CR = outer(ix, ix, function(j, k) 0.8^abs(k - j) * vR[pmin(j, k)]),
       g  = AL_TV * cU + AU_TV)                     # loading of h_k on U
}
#' Mean and sd of h_k under arm `a`.  (Marginal; used for CHECK (c) and the closed form.)
haz_moments_tv <- function(a, alpha.0) {
  M <- .moments_tv(a); k <- seq_len(NV_TV)
  list(mean = alpha.0 + (AA_TV + BETA_TV * (k - 1)) * a + AL_TV * M$b,
       sd   = sqrt(AL_TV^2 * diag(M$CR) + M$g^2 * SIGU_TV^2),
       Sigma = AL_TV^2 * M$CR + outer(M$g, M$g) * SIGU_TV^2, M = M)
}
#' P(h_k <= 0) for each k under arm `a` -- how hard the guard bites.
p_haz_neg_tv <- function(a, alpha.0) {
  H <- haz_moments_tv(a, alpha.0); pnorm(0, H$mean, H$sd)
}
#' UNGUARDED lognormal-MGF survival, E[exp(-sum w_k h_k)].  This is the WRONG answer for the
#' generator whenever the guard bites, and is used only in CHECK (b), at a configuration where
#' P(h_k <= 0) is astronomically small so that guarded and unguarded must agree.
surv_closedform_tv <- function(t, a, alpha.0) {
  if (t <= 0) return(1)
  w <- .wts_tv(t); H <- haz_moments_tv(a, alpha.0)
  exp(-sum(w * H$mean) + drop(w %*% H$Sigma %*% w) / 2)
}
.wts_tv <- function(t) {
  if (t <= 0) return(numeric(NV_TV))
  Kt <- ceiling(t - 1e-12); s <- t - (Kt - 1)
  c(rep(1, Kt - 1), s, rep(0, NV_TV - Kt))[seq_len(NV_TV)]
}

# ---- section (D): forward filter over the U-free Markov chain M -----------------------
#' Precompute, for arm `a` and alpha0, the survival-weighted sub-densities p_K(l | U=u) of
#' M_K on a grid, for K = 1..NV, together with the U nodes/weights.  All t inside interval K
#' reuse p_K, so the whole t-integration costs no extra recursions.
.filter_cache <- new.env(parent = emptyenv())
.filter_tv <- function(a, alpha.0, nL = 3200L, nU = 32L, lo = -16, hi = 14) {
  ky <- paste(a, alpha.0, nL, nU, lo, hi, sep = "|")
  if (!is.null(.filter_cache[[ky]])) return(.filter_cache[[ky]])
  l  <- seq(lo, hi, length.out = nL); dl <- l[2] - l[1]
  wl <- rep(dl, nL); wl[1] <- wl[nL] <- dl / 2          # composite trapezoid
  gu <- gh_rule_tv(nU); u <- SIGU_TV * gu$x             # U ~ N(0, SIGU^2)
  M  <- .moments_tv(a)
  A  <- alpha.0 + (AA_TV + BETA_TV * (seq_len(NV_TV) - 1)) * a
  # h_k(l, u) = A_k + AL*l + g_k*u  -- an nL x nU matrix per k
  hp <- lapply(seq_len(NV_TV), function(k)
    pmax(outer(A[k] + AL_TV * l, M$g[k] * u, "+"), 0))
  P <- vector("list", NV_TV)
  P[[1]] <- matrix(dnorm(l), nL, nU)                    # M_1 ~ N(0,1), free of u
  for (k in seq_len(NV_TV - 1L)) {
    # M_{k+1} = 0.8 M_k + shift + e,  shift = -a + 0.1*k   (U-free: the point of the split)
    Tk <- outer(l, 0.8 * l + (-a + 0.1 * k), function(x, y) dnorm(x - y))
    P[[k + 1L]] <- Tk %*% (P[[k]] * exp(-hp[[k]]) * wl)
  }
  out <- list(l = l, wl = wl, hp = hp, P = P, uw = gu$w, nL = nL, nU = nU)
  .filter_cache[[ky]] <- out
  out
}

#' S_a(t) and f_a(t) = -dS_a/dt, exact for the GUARDED generator.
#' `flt` is a cached .filter_tv() object; pass one to reuse it across many t.
surv_exact_tvfx <- function(t, a, alpha.0, flt = NULL, want_dens = FALSE,
                            nL = 3200L, nU = 32L) {
  if (t <= 0) return(if (want_dens) c(S = 1, f = NA_real_) else 1)
  if (is.null(flt)) flt <- .filter_tv(a, alpha.0, nL = nL, nU = nU)
  K <- ceiling(t - 1e-12); s <- t - (K - 1)
  h <- flt$hp[[K]]; E <- flt$P[[K]] * exp(-s * h)
  S <- sum(flt$uw * colSums(E * flt$wl))
  if (!want_dens) return(S)
  c(S = S, f = sum(flt$uw * colSums(E * h * flt$wl)))
}

#' Exact E1 target: population Lin-Ying constant additive effect (section (E)).
#' Gauss-Legendre inside each unit interval (S and f are analytic in s there; the only
#' non-smoothness in t is at the integers, where the visit -- and hence w -- changes).
gamma_exact_tvfx <- function(alpha.0, tau = NV_TV, nL = 3200L, nU = 32L, nGL = 32L) {
  gl <- gl_rule_tv(nGL)
  f1 <- .filter_tv(1, alpha.0, nL = nL, nU = nU)
  f0 <- .filter_tv(0, alpha.0, nL = nL, nU = nU)
  num <- den <- 0
  for (K in seq_len(tau)) {
    tt <- (K - 1) + (gl$x + 1) / 2; wq <- gl$w / 2
    for (j in seq_along(tt)) {
      B1 <- surv_exact_tvfx(tt[j], 1, alpha.0, flt = f1, want_dens = TRUE)
      B0 <- surv_exact_tvfx(tt[j], 0, alpha.0, flt = f0, want_dens = TRUE)
      S1 <- B1[["S"]]; S0 <- B0[["S"]]; d1 <- B1[["f"]]; d0 <- B0[["f"]]
      num <- num + wq[j] * (S0 * d1 - S1 * d0) / (S1 + S0)
      den <- den + wq[j] * (S1 * S0) / (S1 + S0)
    }
  }
  num / den
}

#' INDEPENDENT reimplementation of S_a(t): tensor Gauss-Hermite on the Cholesky factor of the
#' full 5-D Cov(h), with no Markov structure, no U/M split and no l-grid.  Shares nothing with
#' .filter_tv() but the constants, so agreement is evidence, not tautology.  Converges only
#' algebraically here because of the kink -- that is the point of (D), and is visible in the
#' refinement in m.
surv_gh_tvfx <- function(t, a, alpha.0, m = 40L) {
  if (t <= 0) return(1)
  w <- .wts_tv(t); d <- ceiling(t - 1e-12)
  H <- haz_moments_tv(a, alpha.0)
  S <- H$Sigma[seq_len(d), seq_len(d), drop = FALSE]
  Lc <- t(chol(S)); r <- gh_rule_tv(m)
  idx <- as.matrix(expand.grid(rep(list(seq_len(m)), d)))
  Z <- matrix(r$x[idx], nrow = nrow(idx), ncol = d)
  W <- exp(rowSums(matrix(log(r$w[idx]), nrow = nrow(idx), ncol = d)))
  Hm <- pmax(sweep(Z %*% t(Lc), 2, H$mean[seq_len(d)], "+"), 0)
  sum(W * exp(-drop(Hm %*% w[seq_len(d)])))
}

#' Empirical Lin-Ying two-sample estimator; algebraically the same solution
#' aalen(Surv(Ts,ev) ~ const(A), n.sim=0)$gamma[1,1] targets for time-fixed A with the baseline
#' profiled out.  NOT claimed as a bit-exact identity -- see the driver, which measures the gap.
ly_two_sample_tv <- function(Ts, ev, A) {
  o <- order(Ts); Ts <- Ts[o]; ev <- ev[o]; A <- A[o]; n <- length(Ts)
  n1 <- sum(A == 1); n0 <- n - n1
  cum1 <- n1 - c(0, cumsum(A == 1)[-n]); cum0 <- n0 - c(0, cumsum(A == 0)[-n])
  e <- ev == 1; Abar <- cum1[e] / (cum1[e] + cum0[e])
  tprev <- c(0, Ts[-n])
  sum(A[e] - Abar) / sum((Ts - tprev) * cum1 * cum0 / pmax(cum1 + cum0, 1))
}

# ================================ driver / verification ================================
if (sys.nframe() == 0L) {
  A0 <- c(0.15, 0.25, 0.40)
  akey <- function(a) sprintf("%.2f", as.numeric(a))     # "0.15","0.25","0.40" -- never format()

  cat("=============== CHECK (c): how hard does the haz>0 guard bite? ===============\n")
  cat("If this is small, the whole framing of (A) is unnecessary. It is not small.\n\n")
  cat(sprintf("%-7s %4s %14s %6s   %s\n", "alpha0", "arm", "max P(h_k<=0)", "at k", "P(h_k<=0), k=1..5"))
  for (a0 in A0) for (a in c(1, 0)) {
    p <- p_haz_neg_tv(a, a0)
    cat(sprintf("%-7s %4d %14.4f %6d   %s\n", akey(a0), a, max(p), which.max(p),
                paste(sprintf("%.4f", p), collapse = " ")))
  }

  cat("\n=============== CHECK (b): guard-free configuration vs the closed form ===============\n")
  cat("Separates 'the Gaussian algebra is right' from 'the guard handling is right'.\n")
  cat("Where P(h_k<=0) is negligible, max(h,0) == h and the filter MUST reproduce the\n")
  cat("elementary lognormal MGF exp(-E[sum w h] + Var[sum w h]/2).  This exercises the FULL\n")
  cat("5-step recursion, the transition kernel, the l-grid and the U rule -- and can fail.\n\n")
  cat(sprintf("%-7s %4s %5s %12s %22s %22s %11s\n",
              "alpha0", "arm", "t", "max P(h<=0)", "filter", "closed form", "rel diff"))
  okb <- TRUE
  for (a0 in c(0.40, 0.25)) for (a in c(0, 1)) {
    pmax_ <- max(p_haz_neg_tv(a, a0))
    for (tt in c(1, 2.5, 5)) {
      sf <- surv_exact_tvfx(tt, a, a0, nL = 3200L, nU = 32L)
      sc <- surv_closedform_tv(tt, a, a0)
      rd <- (sf - sc) / sc
      if (pmax_ < 1e-8) okb <- okb && abs(rd) < 5e-8
      cat(sprintf("%-7s %4d %5.2f %12.2e %22.15f %22.15f %+11.2e%s\n",
                  akey(a0), a, tt, pmax_, sf, sc, rd,
                  if (pmax_ < 1e-8) "" else "   (guard bites: must differ)"))
    }
  }
  cat(sprintf("\n  verdict on the guard-free rows: %s\n", if (okb) "PASS" else "FAIL"))

  cat("\n=============== CHECK (e): l-grid and U-rule convergence ===============\n")
  cat("The l-integrand is KINKED, so refinement is the only honest evidence.\n\n")
  cat(sprintf("%-7s %6s %5s %20s %20s\n", "alpha0", "nL", "nU", "S_1(5)", "gamma"))
  for (a0 in A0) {
    for (nl in c(400L, 800L, 1600L, 3200L)) {
      s1 <- surv_exact_tvfx(5, 1, a0, nL = nl, nU = 32L)
      g  <- gamma_exact_tvfx(a0, nL = nl, nU = 32L, nGL = 24L)
      cat(sprintf("%-7s %6d %5d %20.14f %20.14f\n", akey(a0), nl, 32L, s1, g))
    }
    cat("\n")
  }
  cat("U-rule (Gauss-Hermite over U) refinement at nL=1600, S_1(5):\n")
  for (a0 in A0) {
    v <- sapply(c(16L, 24L, 32L, 48L), function(q) surv_exact_tvfx(5, 1, a0, nL = 1600L, nU = q))
    cat(sprintf("  alpha0=%s  ", akey(a0))); cat(sprintf("%.14f ", v)); cat("\n")
  }
  cat("l-range truncation (nL=3200, nU=32), S_1(5) -- rules out a tail-truncation artefact:\n")
  for (a0 in A0) {
    v <- sapply(list(c(-16, 14), c(-20, 18), c(-24, 22)), function(r)
      surv_exact_tvfx(5, 1, a0, flt = .filter_tv(1, a0, nL = 3200L, nU = 32L, lo = r[1], hi = r[2])))
    cat(sprintf("  alpha0=%s  ", akey(a0))); cat(sprintf("%.14f ", v)); cat("\n")
  }
  cat("Gauss-Legendre refinement in t (nL=3200, nU=32):\n")
  for (a0 in A0) {
    v <- sapply(c(8L, 16L, 24L, 32L), function(q)
      gamma_exact_tvfx(a0, nL = 3200L, nU = 32L, nGL = q))
    cat(sprintf("  alpha0=%s  ", akey(a0))); cat(sprintf("%.13f ", v)); cat("\n")
  }

  cat("\n=============== CHECK (d): independent tensor Gauss-Hermite ===============\n")
  cat("Same target, disjoint implementation (full 5-D Cov + Cholesky, no Markov split).\n")
  cat("Kink => algebraic, not spectral, convergence in m -- shown, not hidden.\n\n")
  cat(sprintf("%-7s %4s %20s %14s %14s %14s\n", "alpha0", "arm", "filter S(5)",
              "GH m=8", "GH m=12", "GH m=18"))
  for (a0 in A0) for (a in c(1, 0)) {
    sf <- surv_exact_tvfx(5, a, a0, nL = 3200L, nU = 32L)
    gh <- sapply(c(8L, 12L, 18L), function(m) surv_gh_tvfx(5, a, a0, m = m))
    cat(sprintf("%-7s %4d %20.14f %14.9f %14.9f %14.9f\n", akey(a0), a, sf, gh[1], gh[2], gh[3]))
  }

  cat("\n=============== CHECK: f = -dS/dt against a central difference ===============\n")
  cat("f is differentiated under the integral; verify it against S numerically (can fail).\n\n")
  for (a0 in A0) for (a in c(1, 0)) {
    fl <- .filter_tv(a, a0, nL = 3200L, nU = 32L)
    tt <- 3.4; hstep <- 1e-5
    fa <- surv_exact_tvfx(tt, a, a0, flt = fl, want_dens = TRUE)[["f"]]
    fd <- -(surv_exact_tvfx(tt + hstep, a, a0, flt = fl) -
            surv_exact_tvfx(tt - hstep, a, a0, flt = fl)) / (2 * hstep)
    cat(sprintf("  alpha0=%s a=%d t=%.1f  analytic %.14f  central-diff %.14f  diff %+.2e\n",
                akey(a0), a, tt, fa, fd, fa - fd))
  }

  cat("\n=============== EXACT gamma vs the frozen single draw ===============\n")
  gx <- vapply(A0, gamma_exact_tvfx, 0, nL = 3200L, nU = 32L, nGL = 32L)
  suppressMessages({library(timereg); library(survival)})
  # The study's OWN true_gamma(), reproduced verbatim from dev/robustness_tvfx.R:31-77.
  NV <- NV_TV; AA <- AA_TV; BETA <- BETA_TV; AL <- AL_TV; AU <- AU_TV
  tv_haz <- function(k, A, L, U, alpha.0)
    alpha.0 + (AA + BETA * (k - 1)) * A[, k] + AL * L[, k] + AU * U
  event_times <- function(U, L, A, alpha.0) {
    n <- length(U); Tt <- rep(NA_real_, n)
    for (k in 1:NV) { u <- runif(n); h <- tv_haz(k, A, L, U, alpha.0); nt <- -log(u) / h
      Tt <- ifelse(is.na(Tt) & nt < 1 & h > 0, k - 1 + nt, Tt) }
    list(Ts = ifelse(is.na(Tt), NV, Tt), D = ifelse(is.na(Tt), 0L, 1L))
  }
  pseudo_rct_tv <- function(alpha.0, n = 3e4, seed = 2024) {
    set.seed(seed)
    arm <- function(av) {
      U <- rnorm(n, 0, 0.1); L <- matrix(nrow = n, ncol = NV); L[, 1] <- rnorm(n, U, 1)
      for (k in 2:NV) L[, k] <- rnorm(n, 0.8 * L[, k - 1] - av + 0.1 * (k - 1) + U, 1)
      event_times(U, L, matrix(av, n, NV), alpha.0)
    }
    a1 <- arm(1); a0 <- arm(0)
    data.frame(Ts = c(a1$Ts, a0$Ts), ev = c(a1$D, a0$D), A = rep(c(1, 0), each = n))
  }
  cat("  (first: does ly_two_sample_tv reproduce timereg::aalen const(A)? measured, not asserted)\n")
  for (a0 in A0) {
    d <- pseudo_rct_tv(a0, n = 4000, seed = 11)
    g1 <- unname(aalen(Surv(Ts, ev) ~ const(A), data = d, n.sim = 0)$gamma[1, 1])
    g2 <- ly_two_sample_tv(d$Ts, d$ev, d$A)
    cat(sprintf("   alpha0=%s aalen %.12f  hand %.12f  diff %+.2e\n", akey(a0), g1, g2, g2 - g1))
  }
  frozen <- vapply(A0, function(a0) {
    d <- pseudo_rct_tv(a0, n = 3e4, seed = 2024)
    unname(aalen(Surv(Ts, ev) ~ const(A), data = d, n.sim = 0)$gamma[1, 1])
  }, 0)
  cat(sprintf("\n%-8s %18s %18s %14s\n", "alpha0", "EXACT", "frozen(1 draw)", "displacement"))
  for (i in seq_along(A0))
    cat(sprintf("%-8s %+18.10f %+18.10f %+14.6f\n", akey(A0[i]), gx[i], frozen[i], frozen[i] - gx[i]))

  stopifnot(length(gx) == length(A0), all(is.finite(gx)), all(gx != 0))
  saveRDS(setNames(as.list(gx), akey(A0)), "dev/truth_tvfx_exact.rds")
  cat("\nwrote dev/truth_tvfx_exact.rds  (names: ", paste(akey(A0), collapse = ", "), ")\n", sep = "")

  if ("--mc" %in% commandArgs(TRUE)) {
    cat("\n=============== CHECK (a)/(f): multi-seed Monte Carlo, the study's OWN path ===============\n")
    NSEED <- as.integer(Sys.getenv("NSEED", "600"))
    for (i in seq_along(A0)) {
      v <- vapply(seq_len(NSEED), function(s) {
        d <- pseudo_rct_tv(A0[i], n = 3e4, seed = 800000 + 1000 * i + s)
        ly_two_sample_tv(d$Ts, d$ev, d$A)
      }, 0)
      mu <- mean(v); se <- sd(v) / sqrt(NSEED)
      cat(sprintf("alpha0=%s | MC mean %+.8f (MCSE %.2e, nseed=%d, 1-draw sd %.6f) | exact %+.8f | z = %+.2f\n",
                  akey(A0[i]), mu, se, NSEED, sd(v), gx[i], (mu - gx[i]) / se))
      cat(sprintf("          frozen seed-2024 draw %+.8f => displacement %+.6f (%.2f MC SD)\n",
                  frozen[i], frozen[i] - gx[i], abs(frozen[i] - gx[i]) / sd(v)))
    }
  }
}
