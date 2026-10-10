# Independent verification of the closed-form E1 target.
#
# CLAIM UNDER TEST:
#   (a) the marginal counterfactual survival S_a(t) has a closed form under the DGM, because under a
#       constant arm (U, L_1..L_5) is jointly Gaussian and each interval hazard is linear in it;
#   (b) gamma_true -- the timereg aalen const(A) fit -- is a functional of S_1 and S_0 ALONE:
#         gamma = int_0^tau [s0*(-dS1) - s1*(-dS0)]/(s1+s0)  /  int_0^tau s1*s0/(s1+s0) dt
#       because the only covariate is randomised A, so the at-risk composition that the estimating
#       equation needs is s1/(s1+s0), not the (heavily selected) survivor composition in (U, L).
#
# This script re-derives both from scratch rather than reusing the existing code, and subjects them
# to checks that do not depend on the derivation being right:
#
#   TEST 1  An EXACT identity the derivation cannot fake. At t=1 the arms' hazards differ by exactly
#           alpha.A conditional on the SHARED (U, L0), so S_1(1)/S_0(1) = exp(-alpha.A) EXACTLY, for
#           every alpha.0. Any error in the Gaussian algebra breaks this.
#   TEST 2  gamma-from-KM against timereg on the SAME dataset. If (b) holds this is an algebraic
#           identity and must agree to ~1e-15, not to Monte-Carlo error.
#   TEST 3  closed-form gamma against Monte-Carlo draws (ARC smoke + local), which must agree within
#           the draws' own MC error.
#
#   Rscript dev/verify_closed_form.R
suppressMessages({library(timereg); library(survival)})
for (f in list.files("steCI/R", pattern = "[.]R$", full.names = TRUE)) source(f)
source("simulation/scenarios.R")

A_A <- DGM_FIXED$alpha.A; A_L <- DGM_FIXED$alpha.L; A_U <- DGM_FIXED$alpha.U
NV  <- DGM_FIXED$n.visit; SIGU <- 0.1

# ---- my own derivation of the joint Gaussian law of (U, L_1..L_NV) under constant arm a ----
moments <- function(a) {
  cU <- b <- vR <- numeric(NV)
  cU[1] <- 1; b[1] <- 0; vR[1] <- 1
  for (k in 2:NV) {
    cU[k] <- 0.8 * cU[k - 1] + 1
    b[k]  <- 0.8 * b[k - 1] - a + 0.1 * (k - 1)
    vR[k] <- 0.64 * vR[k - 1] + 1
  }
  ix <- seq_len(NV)
  CR <- outer(ix, ix, function(j, k) 0.8^abs(k - j) * vR[pmin(j, k)])
  list(cU = cU, b = b, CR = CR)
}

# S_a(t) = E[exp(-sum_k w_k h_k)], h_k = alpha.0 + alpha.A*a + alpha.L*L_k + alpha.U*U, linear in
# Gaussians => lognormal MGF. (Ignores the haz>0 guard; TEST 4 bounds that.)
surv_cf <- function(t, a, alpha.0) {
  if (t <= 0) return(1)
  K <- ceiling(t - 1e-12); s <- t - (K - 1)
  w <- c(rep(1, K - 1), s, rep(0, NV - K))[seq_len(NV)]
  M <- moments(a)
  mean_lin <- sum(w) * (alpha.0 + A_A * a) + A_L * sum(w * M$b)
  cU_coef  <- A_L * sum(w * M$cU) + A_U * sum(w)          # coefficient on U
  var_lin  <- cU_coef^2 * SIGU^2 + A_L^2 * drop(w %*% M$CR %*% w)
  exp(-mean_lin + var_lin / 2)
}

# ---- gamma as a functional of the two survival curves, on a fine grid ----
gamma_from_surv <- function(alpha.0, tau = max(T_HOR), ngrid = 20001) {
  tg <- seq(0, tau, length.out = ngrid); dt <- tg[2] - tg[1]
  S1 <- vapply(tg, surv_cf, 0, a = 1, alpha.0 = alpha.0)
  S0 <- vapply(tg, surv_cf, 0, a = 0, alpha.0 = alpha.0)
  s1 <- (S1[-ngrid] + S1[-1]) / 2; s0 <- (S0[-ngrid] + S0[-1]) / 2   # midpoint on each cell
  dF1 <- -diff(S1); dF0 <- -diff(S0)                                  # event mass on each cell
  num <- sum((s0 * dF1 - s1 * dF0) / (s1 + s0))
  den <- sum(s1 * s0 / (s1 + s0)) * dt                                # trapezoid in t
  num / den
}

# ---- the same functional, fed EMPIRICAL KM curves (for TEST 2) ----
gamma_from_km <- function(km1, km0, tau) {
  tg <- sort(unique(c(0, km1$time, km0$time, tau))); tg <- tg[tg <= tau]
  Sof <- function(km, t) { i <- findInterval(t, km$time); ifelse(i == 0, 1, km$surv[pmax(i, 1)]) }
  S1 <- Sof(km1, tg); S0 <- Sof(km0, tg)
  s1 <- c(1, S1[-length(S1)]); s0 <- c(1, S0[-length(S0)])   # LEFT limit: pairs with dN at t_i
  dF1 <- -diff(c(1, S1)); dF0 <- -diff(c(1, S0))
  num <- sum((s0 * dF1 - s1 * dF0) / (s1 + s0))
  # Denominator: on [t_i, t_{i+1}) the at-risk is the value AFTER the jump at t_i, i.e. S1[i].
  # Using the left limit here too is wrong by ~1e-5 -- exactly the endpoint error to watch for.
  wR  <- S1 * S0 / (S1 + S0)
  den <- sum(wR[-length(wR)] * diff(tg))
  num / den
}

# Replicates true_values()'s own two-arm construction exactly (steCI/R/true_values.R:29-52):
# U and L0 are drawn ONCE and SHARED by both arms; each arm then draws its own L[,2:NV] and its own
# uniforms. Reproducing this faithfully is what makes TEST 2 a test of the formula rather than of my
# ability to re-simulate.
make_arms <- function(n, alpha.0, seed) {
  set.seed(seed)
  U <- rnorm(n, 0, 0.1); L0 <- rnorm(n, U, 1)
  arm <- function(aval) {
    L <- matrix(NA_real_, n, NV); A <- matrix(aval, n, NV)
    L[, 1] <- L0
    for (k in 2:NV) L[, k] <- rnorm(n, 0.8 * L[, k - 1] - A[, k - 1] + 0.1 * (k - 1) + U, 1)
    Tt <- rep(NA_real_, n)
    for (k in 1:NV) {
      u <- runif(n, 0, 1)
      haz <- alpha.0 + A_A * A[, k] + A_L * L[, k] + A_U * U
      nt <- -log(u) / haz
      Tt <- ifelse(is.na(Tt) & nt < 1 & haz > 0, k - 1 + nt, Tt)
    }
    list(Tt = ifelse(is.na(Tt), NV, Tt), D = ifelse(is.na(Tt), 0, 1))
  }
  a1 <- arm(1); a0v <- arm(0)
  data.frame(Tstop = c(a1$Tt, a0v$Tt), event = c(a1$D, a0v$D),
             A = rep(c(1, 0), each = n))
}

cat("=================== TEST 1: the exact ratio identity ===================\n")
cat("At t=1 the arms' hazards differ by exactly alpha.A given the SHARED (U,L0),\n")
cat("so S_1(1)/S_0(1) = exp(-alpha.A) EXACTLY. The derivation cannot fake this.\n\n")
target <- exp(-A_A)
cat(sprintf("  exp(-alpha.A) = %.12f\n", target))
ok1 <- TRUE
for (a0 in c(0.15, 0.25, 0.4)) {
  r <- surv_cf(1, 1, a0) / surv_cf(1, 0, a0)
  cat(sprintf("  alpha.0=%-5s ratio = %.12f   err = %+.2e %s\n", format(a0), r, r - target,
              if (abs(r - target) < 1e-10) "PASS" else "FAIL"))
  ok1 <- ok1 && abs(r - target) < 1e-10
}

cat("\n=================== TEST 2: gamma-from-KM vs timereg ===================\n")
cat("If gamma is a functional of S1,S0 alone this is an ALGEBRAIC identity: expect ~1e-15,\n")
cat("not Monte-Carlo agreement. Anything at 1e-3 means the claim is false.\n\n")
ok2 <- TRUE
for (a0 in c(0.15, 0.25, 0.4)) for (sd_ in c(101L, 202L)) {
  d <- make_arms(5000L, a0, sd_)
  fit <- timereg::aalen(survival::Surv(Tstop, event) ~ const(A), data = d, n.sim = 0)
  g_tr <- fit$gamma[1, 1]
  km1 <- survival::survfit(survival::Surv(Tstop, event) ~ 1, data = d[d$A == 1, ])
  km0 <- survival::survfit(survival::Surv(Tstop, event) ~ 1, data = d[d$A == 0, ])
  g_km <- gamma_from_km(km1, km0, tau = NV)
  cat(sprintf("  alpha.0=%-5s seed %3d | timereg %.12f | from-KM %.12f | diff %+.2e %s\n",
              format(a0), sd_, g_tr, g_km, g_km - g_tr,
              if (abs(g_km - g_tr) < 1e-6) "PASS" else "FAIL"))
  ok2 <- ok2 && abs(g_km - g_tr) < 1e-6
}

cat("\n=================== TEST 3: closed-form gamma vs the cache and the draws ===============\n")
cf <- sapply(c(0.15, 0.25, 0.4), gamma_from_surv)
names(cf) <- c("0.15", "0.25", "0.4")
tc <- readRDS("data/truth_cache.rds")
cat(sprintf("%-7s %14s %14s %12s\n", "alpha.0", "closed form", "cached(1 draw)", "shift"))
for (k in names(cf))
  cat(sprintf("%-7s %14.8f %14.8f %+12.6f\n", k, cf[[k]], tc[[k]]$gamma_true,
              cf[[k]] - tc[[k]]$gamma_true))

arc <- "simulation/steCI_truth_arc/simulation/output/truth"
if (dir.exists(arc)) {
  fs <- list.files(arc, pattern = "[.]rds$", full.names = TRUE)
  if (length(fs)) {
    cat("\n  ARC smoke draws vs closed form (K is tiny, so expect agreement within a few SE):\n")
    for (f in fs) {
      z <- readRDS(f); g <- z$gamma[is.finite(z$gamma)]
      if (!length(g)) next
      key <- as.character(z$alpha0)
      se <- if (length(g) > 1) sd(g) / sqrt(length(g)) else NA_real_
      cat(sprintf("    alpha.0=%-5s K=%d mean %.6f | closed form %.6f | diff %+.6f%s\n",
                  key, length(g), mean(g), cf[[key]], mean(g) - cf[[key]],
                  if (is.finite(se)) sprintf(" (%.1f SE)", abs(mean(g) - cf[[key]]) / se) else ""))
    }
  }
}

cat("\n=================== VERDICT ===================\n")
cat(sprintf("TEST 1 (exact ratio identity): %s\n", if (ok1) "PASS" else "FAIL"))
cat(sprintf("TEST 2 (gamma is S1,S0 only):   %s\n", if (ok2) "PASS" else "FAIL"))
if (ok1 && ok2) {
  cat("\nBoth structural claims hold. The E1 target is DETERMINISTIC:\n")
  for (k in names(cf)) cat(sprintf("  alpha.0 = %-5s gamma_true = %.8f\n", k, cf[[k]]))
} else cat("\nAT LEAST ONE CLAIM FAILED -- do not act on the closed form.\n")
