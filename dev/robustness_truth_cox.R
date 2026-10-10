# EXACT-TO-TOLERANCE truths for robustness DGM (i): the COX-MISSPECIFIED generator
# (dev/robustness_misspec.R, Web Table S7).  Replaces the single-frozen-draw
# pseudo-RCT truth `true_gamma()` with a deterministic quadrature evaluation.
#
# ---------------------------------------------------------------------------
# (A) NO CLOSED FORM EXISTS.  Under a sustained arm a, the confounder recursion
#     L_1 = U + e_1,  L_k = 0.8 L_{k-1} - a + 0.1(k-1) + U + e_k,  e_k ~ N(0,1)
#     is linear with Gaussian innovations, so (U, L_1..L_5) is jointly Gaussian and
#     the Cox log-hazard  eta_k = BA*a + BL*L_k + BU*U  is jointly Gaussian too.
#     But the hazard is MULTIPLICATIVE, h_k = alpha0 * exp(eta_k), so
#         S_a(t) = E[ exp(-alpha0 * sum_k w_k exp(eta_k)) ]
#     is the Laplace transform, at argument alpha0, of a SUM OF CORRELATED
#     LOGNORMALS.  The lognormal has no closed-form Laplace transform (its moment
#     problem is indeterminate and E[exp(-c X)] for X lognormal is a genuinely
#     non-elementary function of c), and no linear combination of the eta_k
#     collapses the sum for t > 1.  So the exp(-mu + v/2) lognormal-MGF identity
#     used for the PRIMARY (additive, `alpha0 + aA*A + aL*L + aU*U`) DGM does NOT
#     carry over.  CONFIRMED: no closed form.
#
#     (For t <= 1 only one term survives and S_a(t) = E[exp(-alpha0 t e^{eta_1})]
#     is a 1-D integral -- still not elementary, but computable to machine
#     precision by 1-D Gauss-Hermite.  Used below as an independent check.)
#
# (B) WHAT IS AVAILABLE: a K-dimensional Gaussian expectation, K = ceiling(t) <= 5,
#     with KNOWN mean and covariance.  Evaluated by TENSOR-PRODUCT GAUSS-HERMITE on
#     Cholesky-transformed coordinates.  Justification for tensor GH over sparse
#     grid / QMC: (i) the integrand z |-> exp(-alpha0 sum w_k exp(mu_k + (Lz)_k)) is
#     entire in each coordinate and, crucially, has a very SMALL effective scale --
#     sd(eta_k) <= 0.09 because BL = BU = 0.05 and sd(U) = 0.1 -- so it is nearly
#     affine over the support and GH converges super-geometrically; (ii) K <= 5, so
#     m^K is affordable (m=14 -> 5.4e5 nodes); (iii) tensor GH is deterministic and
#     positive-weighted, giving a monotone, auditable refinement sequence rather
#     than a randomised error estimate.  Convergence is demonstrated below by
#     refining m and reporting successive differences.
#
# (C) F2 SURVIVES MISSPECIFICATION.  gamma_true is DEFINED as the population
#     solution of the Lin-Ying constant-additive-effect estimating equation fitted
#     to the pseudo-RCT, i.e. `aalen(Surv(Ts,ev) ~ const(A))`.  Profiling out the
#     nonparametric baseline gives
#         0 = int sum_i (A_i - Abar(t)) [ dN_i(t) - Y_i(t) A_i gamma dt ],
#         Abar(t) = sum_i Y_i(t) A_i / sum_i Y_i(t).
#     With EQUAL allocation, A randomised, time-fixed, independent of (U,L), and
#     censoring administrative at tau=5 in BOTH arms, the population limits are
#         Abar(t) -> S1/(S1+S0),   E[dN | A=a] = f_a dt,   E[Y | A=a] = S_a,
#     so (the 1/2 allocation factors cancelling)
#         gamma = int [S0 f1 - S1 f0]/(S1+S0) dt  /  int S1 S0/(S1+S0) dt.
#     Every ingredient of this argument is a property of the ESTIMATOR and of the
#     pseudo-RCT DESIGN.  Nothing in it refers to how T was generated: the only
#     covariate is A, so the estimating equation never needs E[covariate | at risk]
#     in (U,L).  Misspecification of the outcome model is exactly what makes gamma
#     a "working-model summary" rather than a structural parameter -- it does not
#     touch the algebra.  VERIFIED: F2 holds verbatim for S7.
#
# (D) NO haz>0 GUARD IS NEEDED, AND NONE EXISTS.  cox_haz returns
#     alpha.0 * exp(...) > 0 always.  The primary DGM needed a truncation guard
#     because its additive hazard could go negative; here the exponential link
#     makes positivity automatic.  So, unlike the primary DGM's closed form, this
#     evaluation is an IDENTITY for the generator, not an approximation --
#     the only error is quadrature error, which is quantified below.
#
#   Rscript dev/robustness_truth_cox.R
# ---------------------------------------------------------------------------

NV_COX <- 5L
BA_COX <- -0.3; BL_COX <- 0.05; BU_COX <- 0.05
SIGU_COX <- 0.1

# ---- Golub-Welsch quadrature rules (no package dependency) ------------------
gh_rule <- function(m) {   # rule for E[g(X)], X ~ N(0,1); weights sum to 1
  b <- sqrt(seq_len(m - 1) / 2)
  J <- matrix(0, m, m); J[cbind(1:(m-1), 2:m)] <- b; J[cbind(2:m, 1:(m-1))] <- b
  ev <- eigen(J, symmetric = TRUE); o <- order(ev$values)
  list(x = sqrt(2) * ev$values[o], w = ev$vectors[1, o]^2)
}
gl_rule <- function(m) {                       # Gauss-Legendre on [-1,1]
  i <- seq_len(m - 1); b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, m, m); J[cbind(1:(m-1), 2:m)] <- b; J[cbind(2:m, 1:(m-1))] <- b
  ev <- eigen(J, symmetric = TRUE); o <- order(ev$values)
  list(x = ev$values[o], w = 2 * ev$vectors[1, o]^2)
}

# ---- moments of (eta_1..eta_NV) under a sustained arm `a` --------------------
.moments_cox <- function(a) {
  cU <- b <- vR <- numeric(NV_COX); cU[1] <- 1; b[1] <- 0; vR[1] <- 1
  for (k in 2:NV_COX) {
    cU[k] <- 0.8 * cU[k - 1] + 1
    b[k]  <- 0.8 * b[k - 1] - a + 0.1 * (k - 1)
    vR[k] <- 0.64 * vR[k - 1] + 1
  }
  ix <- seq_len(NV_COX)
  CR <- outer(ix, ix, function(j, k) 0.8^abs(k - j) * vR[pmin(j, k)])
  g  <- BU_COX + BL_COX * cU                    # loading of U on eta_k
  list(mu = BA_COX * a + BL_COX * b,
       Sigma = outer(g, g) * SIGU_COX^2 + BL_COX^2 * CR)
}

#' Cached tensor Gauss-Hermite node set for the first `d` etas under arm `a`.
#' Returns exp(eta) at each node (n x d) and the product weights.
.eta_nodes <- function(a, d, m) {
  M <- .moments_cox(a)
  S <- M$Sigma[seq_len(d), seq_len(d), drop = FALSE]
  L <- t(chol(S))                                    # lower triangular
  r <- gh_rule(m)
  idx <- as.matrix(expand.grid(rep(list(seq_len(m)), d)))
  Z <- matrix(r$x[idx], nrow = nrow(idx), ncol = d)
  W <- exp(rowSums(matrix(log(r$w[idx]), nrow = nrow(idx), ncol = d)))
  E <- exp(sweep(Z %*% t(L), 2, M$mu[seq_len(d)], "+"))
  list(E = E, W = W)
}

#' S_a(t) and f_a(t) = -dS_a/dt, by tensor Gauss-Hermite.
#' `dim_force` (optional) integrates over that many etas even if w is shorter --
#' used only as a self-consistency check of the covariance/Cholesky machinery.
surv_exact_cox <- function(t, a, alpha.0, m = 12L, nodes = NULL, dim_force = NULL,
                           want_dens = FALSE) {
  if (t <= 0) return(if (want_dens) c(S = 1, f = NA_real_) else 1)
  K <- ceiling(t - 1e-12); s <- t - (K - 1)
  d <- if (is.null(dim_force)) K else dim_force
  w <- numeric(d); if (K > 1) w[seq_len(K - 1)] <- 1; w[K] <- s
  nd <- if (is.null(nodes)) .eta_nodes(a, d, m) else nodes
  lin <- drop(nd$E %*% w)
  e   <- exp(-alpha.0 * lin)
  S   <- sum(nd$W * e)
  if (!want_dens) return(S)
  c(S = S, f = alpha.0 * sum(nd$W * nd$E[, K] * e))
}

#' Exact E1 target: population Lin-Ying constant additive effect (see (C) above).
#' Integrates interval-by-interval with Gauss-Legendre (the integrand is analytic
#' inside each unit interval and kinked at the integers, where the visit changes).
gamma_exact_cox <- function(alpha.0, tau = NV_COX, m = 12L, nGL = 24L) {
  gl <- gl_rule(nGL); num <- den <- 0
  for (K in seq_len(tau)) {
    nd1 <- .eta_nodes(1, K, m); nd0 <- .eta_nodes(0, K, m)
    tt <- (K - 1) + (gl$x + 1) / 2; wq <- gl$w / 2
    for (j in seq_along(tt)) {
      A1 <- surv_exact_cox(tt[j], 1, alpha.0, nodes = nd1, want_dens = TRUE)
      A0 <- surv_exact_cox(tt[j], 0, alpha.0, nodes = nd0, want_dens = TRUE)
      S1 <- A1[["S"]]; S0 <- A0[["S"]]; f1 <- A1[["f"]]; f0 <- A0[["f"]]
      num <- num + wq[j] * (S0 * f1 - S1 * f0) / (S1 + S0)
      den <- den + wq[j] * (S1 * S0) / (S1 + S0)
    }
  }
  num / den
}

#' Empirical Lin-Ying two-sample estimator, algebraically identical to
#' aalen(Surv(Ts,ev) ~ const(A), n.sim=0)$gamma[1,1] for time-fixed A with the
#' baseline profiled out (verified to <5e-16 against timereg; see driver).  This is
#' the FINITE-SAMPLE counterpart of gamma_exact_cox and makes the sanity checks
#' O(n log n) instead of timereg's O(n * #event times).
ly_two_sample <- function(Ts, ev, A) {
  o <- order(Ts); Ts <- Ts[o]; ev <- ev[o]; A <- A[o]; n <- length(Ts)
  n1 <- sum(A == 1); n0 <- n - n1
  cum1 <- n1 - c(0, cumsum(A == 1)[-n]); cum0 <- n0 - c(0, cumsum(A == 0)[-n])
  e <- ev == 1; Abar <- cum1[e] / (cum1[e] + cum0[e])
  tprev <- c(0, Ts[-n])
  sum(A[e] - Abar) / sum((Ts - tprev) * cum1 * cum0 / pmax(cum1 + cum0, 1))
}

# ============================ driver / verification =========================
if (sys.nframe() == 0L) {
  ALPHAS <- c(0.15, 0.25, 0.40)

  cat("=== (1) Gauss-Hermite convergence in m (nodes per dimension) ===\n")
  cat("S_a(t) at the hardest point t=5 (full 5-D), and gamma_exact_cox.\n\n")
  cat(sprintf("%-6s %-4s %20s %20s %20s\n", "alpha0", "m", "S_1(5)", "S_0(5)", "gamma"))
  for (a0 in ALPHAS) {
    for (m in c(6L, 8L, 10L, 12L)) {
      s1 <- surv_exact_cox(5, 1, a0, m = m); s0 <- surv_exact_cox(5, 0, a0, m = m)
      g  <- gamma_exact_cox(a0, m = m, nGL = 24L)
      cat(sprintf("%-6.2f %-4d %20.16f %20.16f %20.16f\n", a0, m, s1, s0, g))
    }
    cat("\n")
  }

  cat("=== (1b) Gauss-Legendre convergence in nGL (t-quadrature), m=12 ===\n")
  for (a0 in ALPHAS) {
    v <- sapply(c(8L, 16L, 24L, 32L, 48L), function(q) gamma_exact_cox(a0, m = 12L, nGL = q))
    cat(sprintf("alpha0=%.2f  ", a0)); cat(sprintf("%.14f ", v)); cat("\n")
  }

  cat("\n=== (2) VERIFY (a): independent 1-D Gauss-Hermite at t=1 ===\n")
  cat("At t=1 only eta_1 enters, so S_a(1)=E[exp(-alpha0 e^{eta_1})] is 1-D.\n")
  cat("Compare (i) 1-D GH with 200 nodes, (ii) the 5-D machinery forced to dim 5.\n\n")
  r200 <- gh_rule(200L)
  for (a0 in ALPHAS) for (a in c(1, 0)) {
    M <- .moments_cox(a); sd1 <- sqrt(M$Sigma[1, 1]); mu1 <- M$mu[1]
    ref <- sum(r200$w * exp(-a0 * exp(mu1 + sd1 * r200$x)))
    f5  <- surv_exact_cox(1, a, a0, m = 12L, dim_force = 5L)
    f1  <- surv_exact_cox(1, a, a0, m = 12L)
    cat(sprintf("alpha0=%.2f a=%d | 1-D ref %.16f | 5-D forced %.16f (diff %+.2e) | native %.16f (diff %+.2e)\n",
                a0, a, ref, f5, f5 - ref, f1, f1 - ref))
  }

  cat("\n=== (3) exact gamma vs the shipped single-draw truth ===\n")
  suppressMessages({library(timereg); library(survival)})
  NV <- NV_COX; BA <- BA_COX; BL <- BL_COX; BU <- BU_COX
  cox_haz <- function(k, A, L, U, alpha.0) alpha.0 * exp(BA * A[, k] + BL * L[, k] + BU * U)
  event_times <- function(U, L, A, alpha.0) {
    n <- length(U); Tt <- rep(NA_real_, n)
    for (k in 1:NV) { u <- runif(n); h <- cox_haz(k, A, L, U, alpha.0); nt <- -log(u) / h
      Tt <- ifelse(is.na(Tt) & nt < 1, k - 1 + nt, Tt) }
    list(Ts = ifelse(is.na(Tt), NV, Tt), D = ifelse(is.na(Tt), 0L, 1L))
  }
  pseudo_rct <- function(alpha.0, n = 3e4, seed = 2024) {
    set.seed(seed)
    arm <- function(av) {
      U <- rnorm(n, 0, 0.1); L <- matrix(nrow = n, ncol = NV); L[, 1] <- rnorm(n, U, 1)
      for (k in 2:NV) L[, k] <- rnorm(n, 0.8 * L[, k - 1] - av + 0.1 * (k - 1) + U, 1)
      event_times(U, L, matrix(av, n, NV), alpha.0)
    }
    a1 <- arm(1); a0 <- arm(0)
    data.frame(Ts = c(a1$Ts, a0$Ts), ev = c(a1$D, a0$D), A = rep(c(1, 0), each = n))
  }
  true_gamma_fast <- function(alpha.0, n = 3e4, seed = 2024) {
    d <- pseudo_rct(alpha.0, n, seed); ly_two_sample(d$Ts, d$ev, d$A)
  }

  cat("(3a) ly_two_sample == timereg::aalen const(A)?\n")
  for (a0 in ALPHAS) {
    d <- pseudo_rct(a0, n = 4000, seed = 11)
    g1 <- unname(aalen(Surv(Ts, ev) ~ const(A), data = d, n.sim = 0)$gamma[1, 1])
    cat(sprintf("   alpha0=%.2f aalen %.12f  hand %.12f  diff %+.2e\n",
                a0, g1, ly_two_sample(d$Ts, d$ev, d$A), ly_two_sample(d$Ts, d$ev, d$A) - g1))
  }

  exact <- sapply(ALPHAS, gamma_exact_cox, m = 12L, nGL = 32L)
  ship  <- sapply(ALPHAS, true_gamma_fast)   # seed 2024, n=3e4/arm: THE SHIPPED TRUTH
  cat(sprintf("\n%-8s %18s %18s %14s\n", "alpha0", "EXACT", "shipped(1 draw)", "shift"))
  for (i in seq_along(ALPHAS))
    cat(sprintf("%-8.2f %+18.10f %+18.10f %+14.6f\n", ALPHAS[i], exact[i], ship[i], ship[i] - exact[i]))

  cat("\n=== (4) VERIFY (b): 400 independent pseudo-RCTs, n=3e4/arm ===\n")
  for (i in seq_along(ALPHAS)) {
    v <- vapply(1:400, function(s) true_gamma_fast(ALPHAS[i], seed = 700000 + 1000 * i + s), 0)
    mu <- mean(v); se <- sd(v) / sqrt(length(v))
    cat(sprintf("alpha0=%.2f | MC mean %+.8f (se %.8f; single-draw sd %.6f) | exact %+.8f | z = %+.2f\n",
                ALPHAS[i], mu, se, sd(v), exact[i], (mu - exact[i]) / se))
  }

  saveRDS(setNames(as.list(exact), as.character(ALPHAS)), "dev/truth_misspec_exact.rds")
  cat("\nwrote dev/truth_misspec_exact.rds\n")
}
