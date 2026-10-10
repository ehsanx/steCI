# Bucket C ATTRIBUTION ablation: is the ~0.005 asymptotic bias omitted-covariate/frailty
# misspecification (baseline-only adjustment) rather than origin-pooling?
#
# The closed form (dev/plim_closed_form.R) found, in the strong-confounding low-event cell, an
# asymptotic bias plim-truth ~ +0.007 at tau=5 of which only ~0.0015 is the pooling component
# (single-origin comparator). It ATTRIBUTED the remaining ~0.005 to the fitted model adjusting on
# BASELINE A and L only, while the DGM hazard depends on TIME-VARYING L_k and an unmeasured frailty
# U. This script CONFIRMS that attribution by ablating the DGM: if the ~0.005 is misspecification,
# turning off the outcome's dependence on L and U must COLLAPSE plim-truth to the pooling floor.
#
# Arms (cell 1: alpha.0=0.15, gamma.L=3.0, gamma.0=0):
#   baseline : alpha.U=0.015, alpha.L=0.015  -- the real DGM; must reproduce plim-truth ~ +0.007.
#   frailtyoff: alpha.U=0,    alpha.L=0.015  -- outcome no longer depends on U (still on L).
#   bothoff  : alpha.U=0,     alpha.L=0      -- outcome depends on NEITHER, so the baseline-A/L
#              Aalen fit is correctly specified: plim-truth must collapse to ~the pooling floor.
#
# Truths: baseline from dev/truth_cache_exact.rds; bothoff is EXACT closed form (constant hazards
# h_a = alpha.0 + alpha.A*a, so S_a(t)=exp(-t h_a), mrd(t)=S_0-S_1); frailtyoff is ALSO exact via
# surv_exact(alpha.U=0) from dev/truth_closed_form.R -- all three arms exact (an earlier version
# scored this arm against a Monte-Carlo true_values() draw). NB in all arms L still evolves (with U) and
# treatment still depends on L via gamma.L -- only the OUTCOME's L/U coefficients are ablated.
#
#   Rscript dev/plim_ablation.R
# Env: PLIM_AB_BATCHES (40), PLIM_AB_BATCH_N (20000), PLIM_AB_GRID (4000); PLIM_AB_OUT (output path,
# default dev/plim_ablation.rds).
suppressMessages({library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang)})
ROOT <- "."
for (f in list.files(file.path(ROOT, "steCI", "R"), pattern = "[.]R$", full.names = TRUE)) source(f)
source(file.path(ROOT, "simulation", "scenarios.R"))
source(file.path(ROOT, "dev", "truth_closed_form.R"))   # surv_exact(): exact alpha.U=0 truth for the frailty-off arm
tv <- readRDS(file.path(ROOT, "dev", "truth_cache_exact.rds"))

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
TMAX <- 5; GRID <- geti("PLIM_AB_GRID", 4000L); DT <- TMAX / GRID; THOR <- T_HOR
stopifnot(GRID %% TMAX == 0L)
A0 <- 0.15; GL <- 3.0; G0 <- 0.0; AA <- DGM_FIXED$alpha.A   # cell 1

# --- machinery (identical to dev/plim_closed_form.R, Rao-Blackwellised) -----------------------
scadd <- function(len, idx, val) {
  keep <- idx >= 1L & idx <= len & is.finite(val)
  if (!any(keep)) return(numeric(len))
  r <- rowsum(val[keep], idx[keep]); out <- numeric(len); out[as.integer(rownames(r))] <- r[, 1]; out
}
new_acc <- function() list(gram = replicate(6, numeric(GRID + 1L), simplify = FALSE),
                           dlam = replicate(3, numeric(GRID + 1L), simplify = FALSE))
accumulate <- function(acc, rows) {
  a <- rows$a; b <- rows$b; w <- rows$w; X <- cbind(1, rows$A, rows$L)
  jlo <- pmax(1L, floor(a / DT) + 1L); jhi <- pmin(GRID, floor(b / DT)); ok <- jhi >= jlo
  ent <- list(c(1,1), c(1,2), c(1,3), c(2,2), c(2,3), c(3,3))
  for (k in seq_along(ent)) {
    p <- ent[[k]][1]; q <- ent[[k]][2]; v <- w * X[, p] * X[, q]
    acc$gram[[k]] <- acc$gram[[k]] + scadd(GRID + 1L, jlo[ok], v[ok]) - scadd(GRID + 1L, jhi[ok] + 1L, v[ok])
  }
  for (p in 1:3) { v <- w * X[, p] * rows$h * DT
    acc$dlam[[p]] <- acc$dlam[[p]] + scadd(GRID + 1L, jlo[ok], v[ok]) - scadd(GRID + 1L, jhi[ok] + 1L, v[ok]) }
  acc
}
solve_Bstar <- function(acc) {
  gc <- lapply(acc$gram, function(x) cumsum(x)[1:GRID]); dc <- lapply(acc$dlam, function(x) cumsum(x)[1:GRID])
  dB <- matrix(0, GRID, 3)
  for (j in 1:GRID) {
    G <- matrix(c(gc[[1]][j],gc[[2]][j],gc[[3]][j], gc[[2]][j],gc[[4]][j],gc[[5]][j], gc[[3]][j],gc[[5]][j],gc[[6]][j]),3,3)
    dB[j, ] <- tryCatch(solve(G, c(dc[[1]][j], dc[[2]][j], dc[[3]][j])), error = function(e) c(NA,NA,NA))
  }
  dB[!is.finite(dB)] <- 0; apply(dB, 2, cumsum)
}
standardise <- function(B, L0, thor) {
  Jt <- pmin(GRID, floor(thor / DT))
  vapply(seq_along(thor), function(k) { b0<-B[Jt[k],1]; bA<-B[Jt[k],2]; bL<-B[Jt[k],3]
    mean(exp(-(b0 + bL*L0))) - mean(exp(-(b0 + bA + bL*L0))) }, numeric(1))
}
pool_acc <- function(per, idx) { a <- new_acc()
  for (k in 1:6) a$gram[[k]] <- Reduce(`+`, lapply(per[idx], function(z) z$gram[[k]]))
  for (p in 1:3) a$dlam[[p]] <- Reduce(`+`, lapply(per[idx], function(z) z$dlam[[p]])); a }

prep_arm <- function(n, seed, aU, aL) {
  d <- dat_sim(n=n, n.visit=DGM_FIXED$n.visit, gamma.0=G0, gamma.L=GL, alpha.0=A0, alpha.A=AA,
               alpha.L=aL, alpha.U=aU, seed=seed)
  set.seed(seed); U <- stats::rnorm(n, 0, 0.1)
  ds <- ste_prepare(d, n.visit=DGM_FIXED$n.visit)
  hpos <- pmax(A0 + AA*ds$A + aL*ds$L + aU*U[ds$id], 0)      # true intensity per row
  list(rows = data.frame(a=ds$time.new, b=ds$time.stop.new, A=ds$A.baseline, L=ds$L.baseline,
                         w=ds$ipw.stab.cum, trial=ds$trial, h=hpos),
       L0 = d$L[d$visit == 1])
}

# Per-arm truth.
truth_baseline <- tv[["0.15"]]$mrd
truth_bothoff  <- exp(-A0*THOR) - exp(-(A0+AA)*THOR)          # exact: S_0 - S_1, constant hazards
truth_frailtyoff <- {                                        # EXACT closed form at alpha.U=0
  s1 <- vapply(THOR, surv_exact, 0, a=1, alpha.0=A0, alpha.U=0)
  s0 <- vapply(THOR, surv_exact, 0, a=0, alpha.0=A0, alpha.U=0)
  s0 - s1
}

ARMS <- list(
  baseline  = list(aU=0.015, aL=0.015, truth=truth_baseline),
  frailtyoff= list(aU=0,     aL=0.015, truth=truth_frailtyoff),
  bothoff   = list(aU=0,     aL=0,     truth=truth_bothoff))

B_N <- geti("PLIM_AB_BATCH_N", 20000L); B_B <- geti("PLIM_AB_BATCHES", 40L)
cat(sprintf("=== ablation, cell 1 (a0=%.2f gL=%.1f) | %d x n=%d | GRID=%d ===\n", A0, GL, B_B, B_N, GRID))
res <- list(); t0 <- Sys.time()
for (nm in names(ARMS)) {
  arm <- ARMS[[nm]]; L0 <- numeric(0); per <- perO <- vector("list", B_B)
  for (bt in seq_len(B_B)) {
    p <- prep_arm(B_N, seed = 700000L + bt, aU = arm$aU, aL = arm$aL)
    per[[bt]]  <- accumulate(new_acc(), p$rows)
    perO[[bt]] <- accumulate(new_acc(), p$rows[p$rows$trial == 0, , drop = FALSE])   # single-origin
    if (length(L0) < 400000L) L0 <- c(L0, p$L0)
  }
  ab <- seq_len(B_B)
  mP <- standardise(solve_Bstar(pool_acc(per,  ab)), L0, THOR)
  mS <- standardise(solve_Bstar(pool_acc(perO, ab)), L0, THOR)
  jk <- vapply(ab, function(b) standardise(solve_Bstar(pool_acc(per, setdiff(ab, b))), L0, THOR), numeric(length(THOR)))
  se <- sqrt((B_B-1)/B_B * rowSums((jk - rowMeans(jk))^2))
  res[[nm]] <- data.frame(arm=nm, tau=THOR, truth=arm$truth, plim=mP, bias=mP-arm$truth,
                          se=se, pooling=mP-mS, stringsAsFactors=FALSE)
  cat(sprintf("\n--- %s (alpha.U=%.3f, alpha.L=%.3f) ---\n", nm, arm$aU, arm$aL))
  for (k in seq_along(THOR))
    cat(sprintf("  tau=%d  plim-truth=%+.5f +/-%.5f  pooling=%+.5f\n", THOR[k], mP[k]-arm$truth[k], se[k], mP[k]-mS[k]))
}
out <- do.call(rbind, res)
cat("\n=== HEADLINE: does the bias collapse when the outcome no longer depends on L,U? ===\n")
bt5 <- function(a) out$bias[out$arm==a & out$tau==5]; st5 <- function(a) out$se[out$arm==a & out$tau==5]
cat(sprintf("  tau=5 plim-truth:  baseline %+.5f (+/-%.5f) | frailtyoff %+.5f | bothoff %+.5f (+/-%.5f)\n",
            bt5("baseline"), st5("baseline"), bt5("frailtyoff"), bt5("bothoff"), st5("bothoff")))
cat(sprintf("  collapse baseline->bothoff = %+.5f  (if ~= the ~0.005 non-pooling bulk, attribution CONFIRMED)\n",
            bt5("baseline") - bt5("bothoff")))
saveRDS(out, Sys.getenv("PLIM_AB_OUT", file.path(ROOT, "dev", "plim_ablation.rds")))
cat(sprintf("\n[ablation] wrote dev/plim_ablation.rds (%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))
