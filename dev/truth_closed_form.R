# EXACT truths for the simulation DGM -- no Monte Carlo.
#
# Both estimands' targets are available in closed form, verified independently in
# dev/verify_closed_form.R:
#
#   E2: S_a(t) = E[exp(-sum_k w_k h_k)]. Under a CONSTANT arm the L-recursion is linear with
#       Gaussian innovations, so (U, L_1..L_5) is jointly Gaussian; the interval hazard is linear
#       in it and constant on [k-1,k); so this is a lognormal MGF, exp(-mu + sigma^2/2).
#
#   E1: gamma_true is the timereg aalen const(A) fit, i.e. the Lin-Ying two-sample estimating
#       equation. Its population solution is a functional of S_1 and S_0 ALONE:
#
#         gamma = int_0^tau [s0*(-dS1) - s1*(-dS0)]/(s1+s0)  /  int_0^tau s1*s0/(s1+s0) dt
#
#       The estimating equation would normally need E[covariate | at risk], but the only covariate
#       is A, which is randomised and time-fixed, so the at-risk composition it needs is
#       s1/(s1+s0) -- a pure function of the two curves. The survivor composition in (U,L) is
#       heavily selected and NOT recoverable from S_a, but it never enters.
#
# ACCURACY. The `haz > 0` guard in the generator makes this an approximation rather than an
# identity, but a negligible one: max_k P(haz_k <= 0) is 7.3e-4 at alpha.0=0.15 in the treated arm
# and <= 1.8e-11 in every other arm-by-event-rate cell (the largest is the alpha.0=0.15 CONTROL arm,
# 1.777e-11 at k=5; 1.894e-13 is the alpha.0=0.25 TREATED arm, not the global runner-up), and the
# resulting error in S is <= 5.6e-6 -- three orders below both the E1 bias (~0.0057) and the
# replication SD of the Monte-Carlo target it replaces (0.0025-0.0063).
#
#   Rscript dev/truth_closed_form.R          # prints the exact truths and what they move
suppressMessages({library(survival)})
for (f in list.files("steCI/R", pattern = "[.]R$", full.names = TRUE)) source(f)
source("simulation/scenarios.R")

A_A <- DGM_FIXED$alpha.A; A_L <- DGM_FIXED$alpha.L; A_U <- DGM_FIXED$alpha.U
NV  <- DGM_FIXED$n.visit; SIGU <- 0.1

.moments <- function(a) {
  cU <- b <- vR <- numeric(NV); cU[1] <- 1; b[1] <- 0; vR[1] <- 1
  for (k in 2:NV) {
    cU[k] <- 0.8 * cU[k - 1] + 1
    b[k]  <- 0.8 * b[k - 1] - a + 0.1 * (k - 1)
    vR[k] <- 0.64 * vR[k - 1] + 1
  }
  ix <- seq_len(NV)
  list(cU = cU, b = b, CR = outer(ix, ix, function(j, k) 0.8^abs(k - j) * vR[pmin(j, k)]))
}

#' Exact marginal counterfactual survival under sustained arm `a`.
surv_exact <- function(t, a, alpha.0, alpha.U = A_U) {
  if (t <= 0) return(1)
  K <- ceiling(t - 1e-12); s <- t - (K - 1)
  w <- c(rep(1, K - 1), s, rep(0, NV - K))[seq_len(NV)]
  M <- .moments(a)
  mu  <- sum(w) * (alpha.0 + A_A * a) + A_L * sum(w * M$b)
  cUc <- A_L * sum(w * M$cU) + alpha.U * sum(w)
  v   <- cUc^2 * SIGU^2 + A_L^2 * drop(w %*% M$CR %*% w)
  exp(-mu + v / 2)
}

#' Exact E1 target: the population Lin-Ying constant additive effect.
gamma_exact <- function(alpha.0, tau = NV, ngrid = 40001) {
  tg <- seq(0, tau, length.out = ngrid); dt <- tg[2] - tg[1]
  S1 <- vapply(tg, surv_exact, 0, a = 1, alpha.0 = alpha.0)
  S0 <- vapply(tg, surv_exact, 0, a = 0, alpha.0 = alpha.0)
  s1 <- (S1[-ngrid] + S1[-1]) / 2; s0 <- (S0[-ngrid] + S0[-1]) / 2
  sum((s0 * (-diff(S1)) - s1 * (-diff(S0))) / (s1 + s0)) / (sum(s1 * s0 / (s1 + s0)) * dt)
}

#' Drop-in replacement for compute_truths(): same shape as truth_cache.rds, exact values.
compute_truths_exact <- function(alpha0_levels = c(0.15, 0.25, 0.4), t.hor = T_HOR) {
  tv <- lapply(alpha0_levels, function(a0) {
    s1 <- vapply(t.hor, surv_exact, 0, a = 1, alpha.0 = a0)
    s0 <- vapply(t.hor, surv_exact, 0, a = 0, alpha.0 = a0)
    list(t.hor = t.hor, surv1 = s1, surv0 = s0, survdiff = s1 - s0, mrd = s0 - s1,
         gamma_true = gamma_exact(a0), exact = TRUE)
  })
  names(tv) <- as.character(alpha0_levels); tv
}

if (sys.nframe() == 0L) {
  old <- readRDS("data/truth_cache.rds")
  new <- compute_truths_exact()

  cat("=========== E1 target: gamma_true ===========\n")
  cat(sprintf("%-7s %14s %14s %12s %10s\n", "alpha.0", "EXACT", "cached(1 draw)", "shift", "% of bias"))
  for (k in names(new)) {
    d <- old[[k]]$gamma_true - new[[k]]$gamma_true
    cat(sprintf("%-7s %14.8f %14.8f %+12.6f %9.0f%%\n",
                k, new[[k]]$gamma_true, old[[k]]$gamma_true, d, 100 * abs(d) / 0.0057))
  }

  cat("\n=========== E2 targets: MRD(tau) ===========\n")
  cat("The cached E2 truths come from KM on a 5e5 pseudo-population, so they carry MC error too --\n")
  cat("much smaller relative to E2 quantities than E1's, but not zero. Movement:\n\n")
  cat(sprintf("%-7s %5s %13s %13s %12s\n", "alpha.0", "tau", "EXACT mrd", "cached mrd", "shift"))
  mx <- 0
  for (k in names(new)) for (i in seq_along(T_HOR)) {
    d <- old[[k]]$mrd[i] - new[[k]]$mrd[i]; mx <- max(mx, abs(d))
    cat(sprintf("%-7s %5d %13.8f %13.8f %+12.8f\n", k, T_HOR[i], new[[k]]$mrd[i], old[[k]]$mrd[i], d))
  }
  cat(sprintf("\nmax |E2 shift| = %.2e", mx))
  cat(sprintf("   (E2 bias at tau=5 runs ~0.001-0.029, so this is ~%.1f%% of the smallest)\n",
              100 * mx / 0.001))

  saveRDS(new, "dev/truth_cache_exact.rds")
  cat("\nwrote dev/truth_cache_exact.rds\n")
}
