# Robustness DGM (i): OUTCOME-MODEL MISSPECIFICATION.
# Event times are generated from a COX (proportional-hazards) mechanism with the
# SAME primary single-confounder structure (frailty U; confounder L_k driven by
# past treatment; monotone treatment), then analysed with the UNCHANGED
# additive-hazards STE estimator (const-aalen). The outcome model is therefore
# deliberately wrong: E1's constant hazard difference is a working-model summary,
# and its true value is taken from a large pseudo-randomised trial analysed the
# SAME (misspecified) additive way. Self-contained; runs LOCALLY; weighted-Aalen
# fits only (no nonparametric bootstrap). Deliverables per cell: coverage of
# nominal 95% intervals for M1 (model)/M2 (row-robust)/M3 (subject-cluster) SEs,
# the SE-ratios, and the (C5) within-subject cross-covariance share
# Delta/V_cl = 1 - (SE_row/SE_cl)^2.
#   Rscript dev/robustness_misspec.R
suppressMessages({library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang)})
for (f in list.files("steCI/R", pattern = "[.]R$", full.names = TRUE)) source(f)
expit <- function(x) exp(x) / (1 + exp(x))
NV <- 5L                             # visits 0..4, admin. censoring at t = 5
BA <- -0.3; BL <- 0.05; BU <- 0.05   # Cox log-hazard-ratio coefficients (A, L, U)

# --- structure shared with the primary DGM: frailty U, confounder L, monotone A
gen_ULA <- function(n, gamma.0, gamma.L) {
  A <- L <- matrix(nrow = n, ncol = NV)
  U <- rnorm(n, 0, 0.1)
  L[, 1] <- rnorm(n, U, 1); A[, 1] <- rbinom(n, 1, expit(gamma.0 + gamma.L * L[, 1]))
  for (k in 2:NV) {
    L[, k] <- rnorm(n, 0.8 * L[, k - 1] - A[, k - 1] + 0.1 * (k - 1) + U, 1)
    A[, k] <- ifelse(A[, k - 1] == 1, 1, rbinom(n, 1, expit(gamma.0 + gamma.L * L[, k])))
  }
  list(U = U, L = L, A = A)
}
# proportional-hazards (Cox) per-visit hazard -- always positive, so no truncation
cox_haz <- function(k, A, L, U, alpha.0) alpha.0 * exp(BA * A[, k] + BL * L[, k] + BU * U)

event_times <- function(U, L, A, alpha.0) {
  n <- length(U); Tt <- rep(NA_real_, n)
  for (k in 1:NV) { u <- runif(n); h <- cox_haz(k, A, L, U, alpha.0); nt <- -log(u) / h
    Tt <- ifelse(is.na(Tt) & nt < 1, k - 1 + nt, Tt) }
  list(Ts = ifelse(is.na(Tt), NV, Tt), D = ifelse(is.na(Tt), 0L, 1L))
}

to_long <- function(U, L, A, Ts, D) {
  n <- length(Ts); vn <- NV - 1
  colnames(A) <- paste0("A.", 0:vn); colnames(L) <- paste0("L.", 0:vn)
  dat <- data.frame(id = 1:n, T.obs = Ts, D.obs = D, A, L)
  dl <- reshape(dat, varying = c(paste0("A.", 0:vn), paste0("L.", 0:vn)),
                direction = "long", idvar = "id")
  dl <- dl[order(dl$id, dl$time), ]
  dl$time.stop <- pmin(dl$time + 1, dl$T.obs)
  dl <- dl[dl$time < dl$T.obs, ]
  dl$event <- as.integer(dl$time.stop == dl$T.obs & dl$D.obs == 1)
  dl$visit <- ave(rep(1, nrow(dl)), dl$id, FUN = cumsum)
  dl <- dl %>% group_by(id) %>%
    mutate(Alag1 = replace_na(lag(A), 0), L.baseline = first(L)) %>% ungroup()
  as.data.frame(dl)[, c("id","time","time.stop","event","visit","A","L","Alag1","L.baseline")]
}
dgm <- function(n, alpha.0, gamma.0 = -1, gamma.L = 1.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  g <- gen_ULA(n, gamma.0, gamma.L); e <- event_times(g$U, g$L, g$A, alpha.0)
  to_long(g$U, g$L, g$A, e$Ts, e$D)
}
# E1 fit: sequential trials + IPACW (package), then const-aalen row and cluster SEs
fit_e1 <- function(d) {
  ds <- ste_prepare(d)
  f0 <- aalen(Surv(time.new, time.stop.new, event) ~ const(A.baseline) + L.baseline,
              data = ds, n.sim = 0, weights = ds$ipw.stab.cum, robust = 1)
  f1 <- aalen(Surv(time.new, time.stop.new, event) ~ const(A.baseline) + L.baseline,
              data = ds, n.sim = 0, weights = ds$ipw.stab.cum, robust = 1, clusters = ds$id)
  c(gamma = unname(f0$gamma[1, 1]), se_model = sqrt(f0$var.gamma[1, 1]),
    se_row = sqrt(f0$robvar.gamma[1, 1]), se_cl = sqrt(f1$robvar.gamma[1, 1]))
}
# marginal working-model truth from a large pseudo-RCT (analysed the SAME additive way)
true_gamma <- function(alpha.0, n = 3e4, seed = 2024) {
  set.seed(seed)
  arm <- function(av) {
    U <- rnorm(n, 0, 0.1); L <- matrix(nrow = n, ncol = NV); L[, 1] <- rnorm(n, U, 1)
    for (k in 2:NV) L[, k] <- rnorm(n, 0.8 * L[, k - 1] - av + 0.1 * (k - 1) + U, 1)
    event_times(U, L, matrix(av, n, NV), alpha.0)
  }
  a1 <- arm(1); a0 <- arm(0)
  dd <- data.frame(Ts = c(a1$Ts, a0$Ts), ev = c(a1$D, a0$D), A = rep(c(1, 0), each = n))
  unname(aalen(Surv(Ts, ev) ~ const(A), data = dd, n.sim = 0)$gamma[1, 1])
}

grid <- expand.grid(n = c(300L, 1000L, 5000L), alpha.0 = c(0.15, 0.25, 0.40))
grid$ev <- c("low","med","high")[match(grid$alpha.0, c(0.15, 0.25, 0.40))]
REPS <- as.integer(Sys.getenv("REPS", "500")); z <- qnorm(0.975)
truths <- sapply(c(0.15, 0.25, 0.40), true_gamma); names(truths) <- c("0.15","0.25","0.4")
cat("Cox-misspecified DGM; marginal working-model gamma (additive analysis):",
    paste(sprintf("%s=%.4f", names(truths), truths), collapse = "  "), "\n\n")

res <- list()
for (gi in seq_len(nrow(grid))) {
  cell <- grid[gi, ]; tg <- truths[[as.character(cell$alpha.0)]]
  est <- matrix(NA, REPS, 4, dimnames = list(NULL, c("gamma","se_model","se_row","se_cl")))
  for (r in seq_len(REPS)) {
    d <- tryCatch(dgm(cell$n, cell$alpha.0, seed = 1000L * gi + r), error = function(e) NULL)
    if (is.null(d)) next
    v <- tryCatch(fit_e1(d), error = function(e) NULL); if (!is.null(v)) est[r, ] <- v
  }
  e <- est[is.finite(est[, "gamma"]) & is.finite(est[, "se_cl"]) & est[, "se_cl"] > 0, , drop = FALSE]
  cov <- c(model = mean(abs(e[,"gamma"] - tg) <= z * e[,"se_model"]),
           row   = mean(abs(e[,"gamma"] - tg) <= z * e[,"se_row"]),
           cl    = mean(abs(e[,"gamma"] - tg) <= z * e[,"se_cl"]))
  sdg <- sd(e[, "gamma"])
  ser <- c(model = mean(e[,"se_model"]) / sdg, row = mean(e[,"se_row"]) / sdg, cl = mean(e[,"se_cl"]) / sdg)
  delta_frac <- mean(1 - (e[,"se_row"] / e[,"se_cl"])^2)
  cat(sprintf("n=%-4d ev=%-4s | cov M1/M2/M3=%.3f/%.3f/%.3f | SEr M2/M3=%.2f/%.2f | Delta/Vcl=%.3f | nvalid=%d MCSE~%.3f\n",
      cell$n, cell$ev, cov["model"], cov["row"], cov["cl"], ser["row"], ser["cl"],
      delta_frac, nrow(e), sqrt(0.25 / nrow(e))))
  res[[gi]] <- data.frame(n = cell$n, ev = cell$ev, cov_model = cov["model"], cov_row = cov["row"],
                          cov_cl = cov["cl"], ser_row = ser["row"], ser_cl = ser["cl"],
                          delta_frac = delta_frac, nvalid = nrow(e), row.names = NULL)
}
cat("\n=== SUMMARY (Cox-misspecified DGM) ===\n"); print(do.call(rbind, res), digits = 3, row.names = FALSE)
saveRDS(do.call(rbind, res), "simulation/output/results/robustness_misspec.rds")
cat("\nROBUSTNESS_MISSPEC DONE\n")
