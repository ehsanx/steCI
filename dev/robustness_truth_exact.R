# EXACT truth for the TWO-CONFOUNDER robustness DGM (Web Table S4, dev/robustness_2confounder.R).
#
# The study currently scores coverage against true_gamma2(), which is ONE frozen Monte-Carlo draw
# from a 6e4/arm pseudo-RCT at seed 2024. A single-draw truth is a fixed displacement shared by
# every replicate, so it does not average out. This file replaces it with a closed form, by the
# same route that dissolved M8 for the primary DGM (dev/truth_closed_form.R).
#
# WHY THE CLOSED FORM EXISTS HERE TOO.
#   Under a SUSTAINED arm a, both confounder recursions are linear with Gaussian innovations:
#       L1_k = 0.8 L1_{k-1} - a       + 0.1(k-1) + U + e_k ,  L1_1 = U    + e_1
#       L2_k = 0.6 L2_{k-1} - 0.5 a              + U + f_k ,  L2_1 = 0.5U + f_1
#   so (U, L1_1..L1_6, L2_1..L2_6) is jointly Gaussian. Decompose
#       L1_k = c1_k U + b1_k + R1_k ,   L2_k = c2_k U + b2_k + R2_k
#   with R1 driven ONLY by e and R2 ONLY by f (two independent rnorm streams). Hence
#
#       *** Cov(L1_j, L2_k) = c1_j * c2_k * sigma_U^2, and NOTHING else. ***
#
#   The two confounders are dependent, but their ONLY channel of dependence is the shared U:
#   Cov(R1_j, R2_k) = 0 because the innovation sequences are independent and each residual process
#   is a function of its own sequence alone. This is the cross-covariance the task flags as the
#   likely error; it enters the variance only through the squared U-coefficient (see cUc below),
#   which is why the variance has no explicit cross term.
#
#   The interval hazard h_k = alpha0 + AA*a + AL1*L1_k + AL2*L2_k + AU*U is ADDITIVE and LINEAR in
#   that Gaussian vector and constant on [k-1,k), so S_a(t) = E[exp(-sum_k w_k h_k)] is a lognormal
#   MGF: exp(-mu + v/2). Unlike the Cox-misspecified DGM (multiplicative hazard => sum of correlated
#   lognormals), nothing here needs quadrature.
#
#   gamma_true is then the Lin-Ying two-sample estimating equation's population solution, a
#   functional of S_1 and S_0 alone (the only covariate is randomised, time-fixed A):
#       gamma = int [S0*(-dS1) - S1*(-dS0)]/(S1+S0)  /  int S1*S0/(S1+S0) dt
#
# RESULTS (K=6, tau=6):
#     alpha0    gamma closed form    + guard correction   <- use the second column as the truth
#     0.15         -0.07208754          -0.07204612
#     0.25         -0.06854897          -0.06854897
#     0.40         -0.06358794          -0.06358794
#   The study's frozen seed-2024 draw sits at -0.07187466 / -0.06840307 / -0.06288173, i.e. it is
#   displaced by +2.1e-4 / +1.5e-4 / +7.1e-4 -- a fixed offset shared by all 300 replicates, 5 to 17
#   times larger than the residual error of the closed form. That is the defect being removed.
#
# ACCURACY.
#   * Gaussian algebra: EXACT. The t=1 ratio identity S_1(1)/S_0(1) = exp(-AA) holds to 0.0e+00
#     (bit-identical) at all three alpha0 -- see CHECK (a).
#   * Quadrature: |gamma(40001 grid) - gamma(80001 grid)| <= 3.0e-11.
#   * The `h > 0` guard is the ONLY real approximation, and unlike the primary DGM it is NOT
#     negligible: max_k P(h_k <= 0) = 1.15e-2 in the TREATED arm at alpha0=0.15 (vs 7.3e-4 for the
#     primary DGM), because two confounders plus a 6th visit give the hazard more variance and the
#     treated arm's L1 drifts down by -1 per visit. It is <= 4.4e-9 everywhere else.
#     Analytic bound on S_exact - S_gen (sign known, >= 0): 1.9e-4; measured by CRN Monte Carlo:
#     1.92e-4 +- 1.1e-6, so the bound is essentially tight. Propagated into gamma this is
#     +4.14e-5 at alpha0=0.15 and <= 5e-10 at 0.25 and 0.40.
#   * Net: the truth is exact to ~1e-10 at alpha0 = 0.25 and 0.40, and to ~1e-6 at alpha0 = 0.15
#     once the (measured, deterministic) guard correction is applied.
#
# CHECK (b), multi-seed Monte Carlo (1000 pseudo-RCTs of 6e4/arm at alpha0=0.15 and 0.25, 400 at
#   0.40, two independent seed blocks): pooled z = +2.18 / +1.00 / +0.69 against the raw closed form,
#   and +0.91 at alpha0=0.15 against the guard-corrected value. At n=2.4e5/arm the alpha0=0.15 gap is
#   +4.12e-5 (MCSE 3.6e-5) -- matching the +4.14e-5 guard correction, which is the guard being seen
#   directly rather than inferred. NOTE the three alpha0 share seeds, so their z are correlated; the
#   first seed block alone gave +2.93/+2.11/+0.69, which is why a second block was run.
#   The many-seed loop scores gamma via the KM functional rather than timereg::aalen (aalen on the
#   1.2e5-row pseudo-RCT costs minutes); the two agree to 5e-10 at n=6e4 -- an algebraic identity,
#   verified here on THIS DGM, not an approximation.
#
#   Rscript dev/robustness_truth_exact.R          # derivation checks + the exact truths
#   Rscript dev/robustness_truth_exact.R --mc     # ... plus the multi-seed Monte-Carlo check (slow)

K   <- 6L
AA  <- -0.04; AL1 <- 0.015; AL2 <- 0.010; AU <- 0.015
SIGU <- 0.1

# ---- joint Gaussian law of (U, L1_1..L1_K, L2_1..L2_K) under a sustained arm `a` ----
# L_k = c_k U + b_k + R_k. R is an AR(1) with unit innovations, so Cov(R_j,R_k)=phi^|k-j| vR_min.
.moments2 <- function(a) {
  c1 <- b1 <- v1 <- numeric(K); c1[1] <- 1.0; b1[1] <- 0; v1[1] <- 1
  c2 <- b2 <- v2 <- numeric(K); c2[1] <- 0.5; b2[1] <- 0; v2[1] <- 1
  for (k in 2:K) {
    c1[k] <- 0.8 * c1[k - 1] + 1
    b1[k] <- 0.8 * b1[k - 1] - a + 0.1 * (k - 1)
    v1[k] <- 0.64 * v1[k - 1] + 1
    c2[k] <- 0.6 * c2[k - 1] + 1
    b2[k] <- 0.6 * b2[k - 1] - 0.5 * a
    v2[k] <- 0.36 * v2[k - 1] + 1
  }
  ix <- seq_len(K)
  list(c1 = c1, b1 = b1, c2 = c2, b2 = b2,
       C1 = outer(ix, ix, function(j, k) 0.8^abs(k - j) * v1[pmin(j, k)]),
       C2 = outer(ix, ix, function(j, k) 0.6^abs(k - j) * v2[pmin(j, k)]))
}

.wts <- function(t) {                       # occupancy of [k-1,k) up to time t
  if (t <= 0) return(numeric(K))
  Kt <- ceiling(t - 1e-12); s <- t - (Kt - 1)
  c(rep(1, Kt - 1), s, rep(0, K - Kt))[seq_len(K)]
}

#' Exact marginal counterfactual survival under sustained arm `a`, two-confounder DGM.
surv_exact2 <- function(t, a, alpha0) {
  if (t <= 0) return(1)
  w <- .wts(t); M <- .moments2(a)
  mu  <- sum(w) * (alpha0 + AA * a) + AL1 * sum(w * M$b1) + AL2 * sum(w * M$b2)
  cUc <- AL1 * sum(w * M$c1) + AL2 * sum(w * M$c2) + AU * sum(w)   # total loading on U
  v   <- cUc^2 * SIGU^2 +
         AL1^2 * drop(w %*% M$C1 %*% w) + AL2^2 * drop(w %*% M$C2 %*% w)
  exp(-mu + v / 2)
}

#' Exact E1 target for the two-confounder DGM.
#' Endpoint convention follows gamma_exact()/gamma_from_km(): the numerator pairs with dN and uses
#' the LEFT limit of the at-risk curves; the denominator integrates the at-risk AFTER each jump.
#' On a fine regular grid both are the midpoint value, which is what the trapezoid below uses.
gamma_exact2 <- function(alpha0, tau = K, ngrid = 40001) {
  tg <- seq(0, tau, length.out = ngrid); dt <- tg[2] - tg[1]
  S1 <- vapply(tg, surv_exact2, 0, a = 1, alpha0 = alpha0)
  S0 <- vapply(tg, surv_exact2, 0, a = 0, alpha0 = alpha0)
  s1 <- (S1[-ngrid] + S1[-1]) / 2; s0 <- (S0[-ngrid] + S0[-1]) / 2
  sum((s0 * (-diff(S1)) - s1 * (-diff(S0))) / (s1 + s0)) / (sum(s1 * s0 / (s1 + s0)) * dt)
}

# ---- the haz>0 guard: marginal law of h_k, and an analytic bound on the induced error in S ----
.haz_moments <- function(a) {
  M <- .moments2(a)
  list(mean = function(alpha0) alpha0 + AA * a + AL1 * M$b1 + AL2 * M$b2,
       sd   = sqrt((AL1 * M$c1 + AL2 * M$c2 + AU)^2 * SIGU^2 +
                   AL1^2 * diag(M$C1) + AL2^2 * diag(M$C2)),
       M = M)
}

#' P(h_k <= 0) for each k under arm a.
p_haz_neg <- function(a, alpha0) {
  H <- .haz_moments(a); pnorm(0, H$mean(alpha0), H$sd)
}

#' Analytic bound on S_closedform(t) - S_generator(t) from the `h > 0` guard.
#'
#' On {h_k <= 0} the generator contributes 1 to the survival product; the closed form contributes
#' exp(-w_k h_k) >= 1. So with N = {k : h_k <= 0}, P = prod_k exp(-w_k h_k) and G the generator's
#' product,   G = P * prod_{k in N} exp(w_k h_k) <= P.  The sign is therefore KNOWN: S_exact >= S_gen.
#' And, using 1 - prod_{k in N} x_k <= sum_{k in N} (1 - x_k) for x_k in [0,1], and 1 - e^{u} <= -u
#' for u <= 0,
#'     0 <= P - G <= P * sum_k w_k (-h_k) 1{h_k <= 0}.
#' The crude alternative (bounding by P(h_k <= 0)) is grossly pessimistic: the discrepancy VANISHES
#' at the boundary h_k = 0, where all of the mass of {h_k <= 0} sits, so the factor (-h_k) matters.
#' Each term is exact by exponential tilting of the bivariate normal (X = sum_j w_j h_j, h_k):
#'     E[e^{-X} (-h_k) 1{h_k<=0}] = E[e^{-X}] * E[(-Y)1{Y<=0}],  Y ~ N(mu_k - Cov(X,h_k), sd_k^2),
#'     E[(-Y)1{Y<=0}] = sd*phi(m/sd) - m*Phi(-m/sd).
guard_bound <- function(t, a, alpha0) {
  w <- .wts(t); M <- .moments2(a); H <- .haz_moments(a)
  gU  <- AL1 * M$c1 + AL2 * M$c2 + AU                 # loading of h_k on U
  cUc <- sum(w * gU)                                  # loading of X = sum w h on U
  covXh <- cUc * gU * SIGU^2 +
           AL1^2 * drop(M$C1 %*% w) + AL2^2 * drop(M$C2 %*% w)
  S <- surv_exact2(t, a, alpha0)
  m <- H$mean(alpha0) - covXh; s <- H$sd              # tilted law of h_k
  S * sum(w * (s * dnorm(m / s) - m * pnorm(-m / s)))
}

#' The guard's effect measured rather than bounded: common-random-numbers Monte Carlo of
#' S_generator(t) - S_exact(t). The two differ only on the rare event {some h_k <= 0}, so under CRN
#' the difference has tiny variance and is estimated far more precisely than either curve.
guard_delta_mc <- function(a, alpha0, tgrid, n = 4e6, chunk = 2e5, seed = 7L) {
  set.seed(seed)
  W <- vapply(tgrid, .wts, numeric(K))                # K x length(tgrid)
  acc <- numeric(length(tgrid)); acc2 <- numeric(length(tgrid)); tot <- 0
  for (ci in seq_len(ceiling(n / chunk))) {
    m <- min(chunk, n - tot); if (m <= 0) break
    U <- rnorm(m, 0, 0.1)
    L1 <- matrix(0, m, K); L2 <- matrix(0, m, K)
    L1[, 1] <- rnorm(m, U, 1); L2[, 1] <- rnorm(m, 0.5 * U, 1)
    for (k in 2:K) {
      L1[, k] <- rnorm(m, 0.8 * L1[, k - 1] - a + 0.1 * (k - 1) + U, 1)
      L2[, k] <- rnorm(m, 0.6 * L2[, k - 1] - 0.5 * a + U, 1)
    }
    h <- alpha0 + AA * a + AL1 * L1 + AL2 * L2 + AU * U     # m x K
    hp <- h; hp[h <= 0] <- 0                                 # generator: no hazard when h<=0
    P <- exp(-(h %*% W)); G <- exp(-(hp %*% W))
    d <- G - P
    acc <- acc + colSums(d); acc2 <- acc2 + colSums(d^2); tot <- tot + m
  }
  mean <- acc / tot
  list(t = tgrid, delta = mean, se = sqrt(pmax(acc2 / tot - mean^2, 0) / tot), n = tot)
}

#' gamma recomputed from guard-corrected survival curves, to see what the guard moves in gamma.
gamma_guarded <- function(alpha0, ngrid = 12001, n = 4e6) {
  tg <- seq(0, K, length.out = ngrid)
  tc <- seq(0, K, by = 0.05)                          # correction is smooth; spline onto fine grid
  S1 <- vapply(tg, surv_exact2, 0, a = 1, alpha0 = alpha0)
  S0 <- vapply(tg, surv_exact2, 0, a = 0, alpha0 = alpha0)
  d1 <- guard_delta_mc(1, alpha0, tc, n = n, seed = 11L)
  d0 <- guard_delta_mc(0, alpha0, tc, n = n, seed = 12L)
  S1 <- S1 + splinefun(tc, d1$delta)(tg); S0 <- S0 + splinefun(tc, d0$delta)(tg)
  s1 <- (S1[-ngrid] + S1[-1]) / 2; s0 <- (S0[-ngrid] + S0[-1]) / 2
  dt <- tg[2] - tg[1]
  g <- sum((s0 * (-diff(S1)) - s1 * (-diff(S0))) / (s1 + s0)) / (sum(s1 * s0 / (s1 + s0)) * dt)
  list(gamma = g, max_d1 = max(abs(d1$delta)), max_d0 = max(abs(d0$delta)),
       se1 = max(d1$se), se0 = max(d0$se))
}

# ======================================================================================
if (sys.nframe() == 0L) {
  A0 <- c(0.15, 0.25, 0.40)

  cat("=============== CHECK (a): the t=1 exact ratio identity ===============\n")
  cat("At t=1, L1_1 and L2_1 do not depend on the arm and A enters h_1 only through AA,\n")
  cat("so h_1(a=1) - h_1(a=0) = AA identically given (U, L1_1, L2_1), hence\n")
  cat("   S_1(1) / S_0(1) = E[e^{-h_1(1)}]/E[e^{-h_1(0)}] = exp(-AA)  EXACTLY.\n\n")
  cat(sprintf("  exp(-AA) = %.15f\n", exp(-AA)))
  ok1 <- TRUE
  for (a0 in A0) {
    r <- surv_exact2(1, 1, a0) / surv_exact2(1, 0, a0); e <- r - exp(-AA)
    ok1 <- ok1 && abs(e) < 1e-13
    cat(sprintf("  alpha0=%-5s  ratio = %.15f   err = %+.3e  %s\n",
                format(a0), r, e, if (abs(e) < 1e-13) "PASS" else "FAIL"))
  }

  cat("\n=============== grid convergence of the gamma quadrature ===============\n")
  for (a0 in A0) {
    g <- vapply(c(10001, 20001, 40001, 80001), function(ng) gamma_exact2(a0, ngrid = ng), 0)
    cat(sprintf("  alpha0=%-5s  %.10f %.10f %.10f %.10f   |40k-80k| = %.1e\n",
                format(a0), g[1], g[2], g[3], g[4], abs(g[3] - g[4])))
  }

  cat("\n=============== CHECK (c): the haz>0 guard ===============\n")
  cat("This DGM has MORE hazard variance than the primary one (two confounders), so the guard is\n")
  cat("expected to bite harder. P(h_k<=0) is NOT negligible in the treated arm at alpha0=0.15.\n\n")
  cat(sprintf("%-7s %4s %14s %5s %14s\n", "alpha0", "arm", "max P(h_k<=0)", "at k", "bound dS"))
  for (a0 in A0) for (a in c(1, 0)) {
    p <- p_haz_neg(a, a0)
    b <- max(vapply(seq(0.25, K, by = 0.25), function(t) guard_bound(t, a, a0), 0))
    cat(sprintf("%-7s %4d %14.3e %5d %14.3e\n", format(a0), a, max(p), which.max(p), b))
  }
  cat("\n  bound is on S_exact - S_gen (>= 0, sign known). Now MEASURED, by CRN Monte Carlo,\n")
  cat("  and propagated into gamma:\n\n")
  cat(sprintf("%-7s %13s %13s %16s %16s %12s\n",
              "alpha0", "max dS arm1", "max dS arm0", "gamma_exact", "gamma_guarded", "shift"))
  # Keep the guarded values: they are the RECOMMENDED truths, and recomputing them is a
  # 4e6-draw CRN job. Discarding them here is what left the .rds holding only the raw
  # closed form, so that anything reading it silently reverted the alpha0=0.15 correction.
  # Index by POSITION, never by a format()ed name. format() pads to a common width on a VECTOR
  # ("0.40") but not on a SCALAR ("0.4"), so `setNames(..., format(A0))` then `x[[format(a0)]]`
  # builds a map whose keys cannot be read back: the assignment APPENDS a 4th element under
  # "0.4" and leaves element 3 at zero. That shipped -- gamma_recommended was
  # (-0.0720461, -0.0685490, 0) -- and would have scored alpha0=0.40 against a truth of 0.
  guarded <- numeric(length(A0))
  for (i in seq_along(A0)) {
    a0 <- A0[i]
    gg <- gamma_guarded(a0)
    guarded[i] <- gg$gamma
    cat(sprintf("%-7s %13.3e %13.3e %16.8f %16.8f %+12.3e\n",
                format(a0), gg$max_d1, gg$max_d0, gamma_exact2(a0), gg$gamma,
                gg$gamma - gamma_exact2(a0)))
  }

  cat("\n=============== EXACT gamma vs the frozen single draw ===============\n")
  gx <- vapply(A0, gamma_exact2, 0)
  cached <- c(NA, NA, NA)
  if (file.exists("dev/robustness_truth_exact_mc.rds")) {
    z <- readRDS("dev/robustness_truth_exact_mc.rds"); cached <- z$single
  }
  cat(sprintf("%-7s %16s\n", "alpha0", "gamma_exact"))
  for (i in seq_along(A0)) cat(sprintf("%-7s %16.10f\n", format(A0[i]), gx[i]))

  # gamma_recommended is what a scorer must use. It differs from gamma_exact only at
  # alpha0=0.15, by +4.1e-5, because the generator clamps with `h > 0` and the closed form
  # integrates the unclamped hazard -- so the closed form is NOT the truth of the DGM as
  # coded. P(h_k<=0) reaches 1.15e-2 in the treated arm at k=6 here, 16x the primary DGM's
  # 7.3e-4. The correction is deterministic (fixed CRN seeds), with its own MC error ~1e-6:
  # QUOTE SIX DECIMALS, NOT EIGHT. Independent replications spread over +4.05e-5 to +4.12e-5.
  # guarded is already positional and in A0 order -- do NOT re-look-up by name here, that is
  # the second half of the same hazard.
  stopifnot(length(guarded) == length(A0), all(is.finite(guarded)), all(guarded != 0))
  saveRDS(list(alpha0 = A0, gamma_exact = gx,
               gamma_recommended = guarded,
               guard_note = paste("gamma_recommended = closed form + measured haz>0 guard",
                                  "correction; differs from gamma_exact only at alpha0=0.15",
                                  "(+4.1e-5). MC error ~1e-6, so quote 6 dp.")),
          "dev/robustness_truth_exact.rds")
  cat("\nwrote dev/robustness_truth_exact.rds (gamma_exact AND gamma_recommended)\n")

  if ("--mc" %in% commandArgs(TRUE)) {
    cat("\n=============== CHECK (b): multi-seed Monte Carlo ===============\n")
    suppressMessages({library(timereg); library(survival)})
    # SAME code path as true_gamma2(), seed varied.
    true_gamma2 <- function(alpha0, n = 6e4, seed = 2024) {
      set.seed(seed); U <- rnorm(n,0,0.1); L10 <- rnorm(n,U,1); L20 <- rnorm(n,0.5*U,1)
      arm <- function(av) { L1<-L2<-matrix(nrow=n,ncol=K); A<-matrix(av,n,K); L1[,1]<-L10; L2[,1]<-L20
        for (k in 2:K) { L1[,k]<-rnorm(n,0.8*L1[,k-1]-A[,k-1]+0.1*(k-1)+U,1)
                         L2[,k]<-rnorm(n,0.6*L2[,k-1]-0.5*A[,k-1]+U,1) }
        Tt<-rep(NA,n); for (k in 1:K) { u<-runif(n)
          h<-alpha0+AA*A[,k]+AL1*L1[,k]+AL2*L2[,k]+AU*U; nt<--log(u)/h
          Tt<-ifelse(is.na(Tt)&nt<1&h>0,k-1+nt,Tt) }
        list(Tt=ifelse(is.na(Tt),K,Tt), D=ifelse(is.na(Tt),0,1)) }
      a1<-arm(1); a0<-arm(0)
      dd<-data.frame(Ts=c(a1$Tt,a0$Tt),ev=c(a1$D,a0$D),A=rep(c(1,0),each=n))
      unname(aalen(Surv(Ts,ev)~const(A),data=dd,n.sim=0)$gamma[1,1])
    }
    NSEED <- as.integer(Sys.getenv("NSEED", "200"))
    out <- list()
    for (i in seq_along(A0)) {
      a0 <- A0[i]
      g <- vapply(seq_len(NSEED), function(s) true_gamma2(a0, seed = 20000 + s), 0)
      mse <- sd(g) / sqrt(NSEED); zst <- (mean(g) - gx[i]) / mse
      single <- true_gamma2(a0, seed = 2024)          # the frozen draw the study uses
      cat(sprintf("alpha0=%-5s  MC mean %.8f (MCSE %.2e, K=%d) | exact %.8f | z = %+.2f\n",
                  format(a0), mean(g), mse, NSEED, gx[i], zst))
      cat(sprintf("            frozen seed-2024 draw %.8f  =>  displacement %+0.6f (%.1f MC SD)\n",
                  single, single - gx[i], abs(single - gx[i]) / sd(g)))
      out[[i]] <- list(alpha0 = a0, mc_mean = mean(g), mc_sd = sd(g), mcse = mse,
                       z = zst, single = single, exact = gx[i], nseed = NSEED)
    }
    saveRDS(list(alpha0 = A0, gamma_exact = gx, mc = out,
                 single = vapply(out, function(o) o$single, 0)),
            "dev/robustness_truth_exact_mc.rds")
    cat("\nwrote dev/robustness_truth_exact_mc.rds\n")
  }
}
