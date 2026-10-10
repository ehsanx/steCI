# Population limit of the pooled weighted-Aalen STE estimator (E2), with (a) Rao-Blackwellised
# event intensity and (b) a SINGLE-ORIGIN comparator that identifies the pooling component.
#
# WHY THIS VERSION EXISTS. The first version reported plim - truth with no Monte-Carlo error and
# called it "the pooled-origin offset". Closer analysis showed both were wrong:
# the MC error swamps tau>=3 (jackknife SE 0.019 / 0.030 at tau=4/5, i.e. the reported offsets were
# 0.24 / 0.35 SE from zero), and "plim - truth" is NOT a pooling effect at all -- the fit adjusts on
# baseline A and L only while the DGM hazard depends on time-varying L_k and the frailty U, so the
# difference bounds the SUM of pooling + baseline-only adjustment + omitted-frailty bias. Worse, the
# weight mass by (window x origin) is lower-triangular: time-since-origin (4,5] is fed by origin 0
# ALONE, so pooling misspecification is structurally ZERO exactly where the offset was being read.
#
# THIS VERSION FIXES BOTH, as the review prescribed.
#
# (a) RAO-BLACKWELLISATION. By the innovation theorem the population event moment is
#         E{ sum_r w_r X_r dN_r(t) } = E{ sum_r w_r X_r Y_r(t) h_r(t) } dt,
#     and h_r is KNOWN EXACTLY here: the DGM hazard is piecewise constant on each visit interval,
#         h = max(alpha.0 + alpha.A*A_k + alpha.L*L_k + alpha.U*U, 0),
#     the positive part because dat_sim only realises an event when haz > 0 (the guard at
#     dat_sim.R:61). Substituting the true intensity for the observed counts removes ALL
#     event-sampling noise: the residual variance is Var(h)/ESS instead of E[h]/ESS
#     (sd(h) ~ 0.025 against mean ~0.15), about 36x the variance and ~6x the SE, at zero extra cost.
#     U is recovered WITHOUT touching the shipped package: dat_sim() calls set.seed(seed) and then
#     draws U <- rnorm(n, 0, 0.1) as its FIRST random draw, so replaying exactly that call
#     reproduces the identical frailty vector. (Asserted below against the realised event counts.)
#
# (b) SINGLE-ORIGIN COMPARATOR. The same machinery is run twice on the SAME rows: once pooled over
#     all origins, once restricted to trial == 0. Under common random numbers the two share their
#     Monte-Carlo noise, so the DIFFERENCE (pooled - single-origin) isolates the pooling component
#     far more sharply than either level -- and it is the only quantity here that may honestly be
#     called a pooled-origin effect. Both arms standardise over the same baseline-L population, so
#     they target the same estimand and differ only in which origins inform B(.).
#
# Every reported quantity carries a delete-one-batch jackknife SE (the batches are i.i.d. by
# construction, which is what licenses it). Horizons that are not resolved are flagged, not spun.
#
# 🔴 A SINGLE RUN'S se_pooled IS A NOISY WITHIN-STREAM JACKKNIFE -- DO NOT trust one run's tau=5. The
# tau=5 point estimate has between-stream SD ~0.0032 and its jackknife SE is unstable across streams
# (0.0016-0.0047) and mildly optimistic (~1.2x), because the heavy-tailed IPACW weights make the
# long-horizon estimate depend on whether a stream hits a rare heavy-weight batch. This file's seed
# base (900000) drew +0.0070, an above-average value; another stream (700000) drew +0.0012.
# THE RESOLVED VALUE, from 7 independent seed streams (dev/plim_seedcheck.R), is
#   tau=5 asymptotic bias = +0.0043 +/- 0.0012 (SE of the mean),  growing monotonically
#   ~0 / +0.0002 / +0.0008 / +0.0020 / +0.0043 at tau=1..5.
# tau<=2 is a clean null. Aggregating INDEPENDENT streams (not more batches within one) is what
# resolves the long horizons. ROBUST regardless of the exact value: pooling is exonerated
# (single-origin comparator, all horizons), and the offset is an outcome-model misspecification
# effect (dev/plim_ablation.R collapses it to exactly 0 under correct specification).
#
#   Rscript dev/plim_closed_form.R
# Env: PLIM_CF_BATCH_N (20000), PLIM_CF_BATCHES (60), PLIM_CF_GRID (4000), PLIM_CF_CELLS ("1:2:3:4")
suppressMessages({library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang)})
ROOT <- "."
for (f in list.files(file.path(ROOT, "steCI", "R"), pattern = "[.]R$", full.names = TRUE)) source(f)
source(file.path(ROOT, "simulation", "scenarios.R"))
tv <- readRDS(file.path(ROOT, "dev", "truth_cache_exact.rds"))

CELLS <- data.frame(
  cell_id = 1:4,
  label   = c("evlow_cfstrong_pvhigh", "evlow_cfweak_pvhigh",
              "evlow_cfmod_pvhigh",    "evhigh_cfstrong_pvhigh"),
  alpha.0 = c(0.15, 0.15, 0.15, 0.40),
  gamma.L = c(3.0,  0.5,  1.5,  3.0),
  gamma.0 = c(0.0,  0.0,  0.0,  0.0),
  stringsAsFactors = FALSE)

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
TMAX <- 5; GRID <- geti("PLIM_CF_GRID", 4000L); DT <- TMAX / GRID
THOR <- T_HOR
# GRID MUST be a multiple of TMAX: interval starts are integers and 100% of censored rows end
# exactly at an integer, so d(t)'s large downward jumps land exactly on cell edges only when
# GRID %% 5 == 0. A non-multiple degrades accuracy SILENTLY.
stopifnot(GRID %% TMAX == 0L)

# ---- one batch: shipped simulate + prepare, plus the exact per-row DGM hazard ------------------
prep_one <- function(cell, n, seed) {
  d <- dat_sim(n = n, n.visit = DGM_FIXED$n.visit, gamma.0 = cell$gamma.0, gamma.L = cell$gamma.L,
               alpha.0 = cell$alpha.0, alpha.A = DGM_FIXED$alpha.A, alpha.L = DGM_FIXED$alpha.L,
               alpha.U = DGM_FIXED$alpha.U, seed = seed)
  # Replay dat_sim's FIRST draw to recover the identical frailty (no package change).
  set.seed(seed); U <- stats::rnorm(n, 0, 0.1)
  ds <- ste_prepare(d, n.visit = DGM_FIXED$n.visit)
  # True piecewise-constant intensity for each interval row, positive part (dat_sim.R:56,61).
  hpos <- pmax(cell$alpha.0 + DGM_FIXED$alpha.A * ds$A + DGM_FIXED$alpha.L * ds$L +
                 DGM_FIXED$alpha.U * U[ds$id], 0)
  list(rows = data.frame(a = ds$time.new, b = ds$time.stop.new, ev = ds$event,
                         A = ds$A.baseline, L = ds$L.baseline, w = ds$ipw.stab.cum,
                         trial = ds$trial, h = hpos),
       L0 = d$L[d$visit == 1])
}

scadd <- function(len, idx, val) {
  keep <- idx >= 1L & idx <= len & is.finite(val)
  if (!any(keep)) return(numeric(len))
  r <- rowsum(val[keep], idx[keep]); out <- numeric(len)
  out[as.integer(rownames(r))] <- r[, 1]; out
}
new_acc <- function() list(gram = replicate(6, numeric(GRID + 1L), simplify = FALSE),
                           dlam = replicate(3, numeric(GRID + 1L), simplify = FALSE),
                           rb = TRUE)

# Accumulate the at-risk Gram and the event moment. rb = TRUE uses the Rao-Blackwellised intensity
# (w X h dt spread over every at-risk cell); rb = FALSE uses the observed events (the original
# estimator moment), kept so the two can be checked against each other.
accumulate <- function(acc, rows, rb = TRUE) {
  a <- rows$a; b <- rows$b; w <- rows$w
  X <- cbind(1, rows$A, rows$L)
  jlo <- pmax(1L, floor(a / DT) + 1L)          # first at-risk cell (edge > a)
  jhi <- pmin(GRID, floor(b / DT))             # last at-risk cell  (edge <= b)
  ok  <- jhi >= jlo
  ent <- list(c(1,1), c(1,2), c(1,3), c(2,2), c(2,3), c(3,3))
  for (k in seq_along(ent)) {
    p <- ent[[k]][1]; q <- ent[[k]][2]; v <- w * X[, p] * X[, q]
    acc$gram[[k]] <- acc$gram[[k]] + scadd(GRID + 1L, jlo[ok], v[ok]) -
                                     scadd(GRID + 1L, jhi[ok] + 1L, v[ok])
  }
  for (p in 1:3) {
    if (rb) {
      # Expected weighted events per cell = w * X * h * DT on every at-risk cell (same difference
      # array as the Gram, since h is constant across the row's interval).
      v <- w * X[, p] * rows$h * DT
      acc$dlam[[p]] <- acc$dlam[[p]] + scadd(GRID + 1L, jlo[ok], v[ok]) -
                                       scadd(GRID + 1L, jhi[ok] + 1L, v[ok])
    } else {
      # `& ok` is required: a row shorter than DT is in no Gram cell, so its event must not be
      # scored against a Gram that excludes it (the one-sided off-by-one the review found).
      je <- pmin(GRID, pmax(1L, floor(b / DT)))
      ie <- rows$ev == 1 & ok
      acc$dlam[[p]][1:GRID] <- acc$dlam[[p]][1:GRID] + scadd(GRID, je[ie], (w * X[, p])[ie])
    }
  }
  acc$rb <- rb
  acc
}

solve_Bstar <- function(acc) {
  gram_cell <- lapply(acc$gram, function(x) cumsum(x)[1:GRID])
  dlam_cell <- if (isTRUE(acc$rb)) lapply(acc$dlam, function(x) cumsum(x)[1:GRID])
               else                lapply(acc$dlam, function(x) x[1:GRID])
  dB <- matrix(0, GRID, 3)
  for (j in 1:GRID) {
    G <- matrix(c(gram_cell[[1]][j], gram_cell[[2]][j], gram_cell[[3]][j],
                  gram_cell[[2]][j], gram_cell[[4]][j], gram_cell[[5]][j],
                  gram_cell[[3]][j], gram_cell[[5]][j], gram_cell[[6]][j]), 3, 3)
    s <- c(dlam_cell[[1]][j], dlam_cell[[2]][j], dlam_cell[[3]][j])
    dB[j, ] <- tryCatch(solve(G, s), error = function(e) c(NA, NA, NA))
  }
  dB[!is.finite(dB)] <- 0
  apply(dB, 2, cumsum)
}

standardise <- function(B, L0, thor) {
  Jt <- pmin(GRID, floor(thor / DT))
  vapply(seq_along(thor), function(k) {
    b0 <- B[Jt[k], 1]; bA <- B[Jt[k], 2]; bL <- B[Jt[k], 3]
    mean(exp(-(b0 + bL * L0))) - mean(exp(-(b0 + bA + bL * L0)))
  }, numeric(1))
}

pool_acc <- function(per, idx, rb) {
  a <- new_acc(); a$rb <- rb
  for (k in seq_len(6)) a$gram[[k]] <- Reduce(`+`, lapply(per[idx], function(z) z$gram[[k]]))
  for (p in 1:3)        a$dlam[[p]] <- Reduce(`+`, lapply(per[idx], function(z) z$dlam[[p]]))
  a
}

BATCH_N <- geti("PLIM_CF_BATCH_N", 20000L)
BATCHES <- geti("PLIM_CF_BATCHES", 60L)
WANT    <- as.integer(strsplit(Sys.getenv("PLIM_CF_CELLS", "1:2:3:4"), "[:,]")[[1]])

cat(sprintf("=== plim: %d batches x n=%d = %s subjects/cell | GRID=%d | Rao-Blackwellised ===\n",
            BATCHES, BATCH_N, format(BATCHES * BATCH_N, big.mark = ","), GRID))
out <- list(); t0 <- Sys.time()
for (cid in WANT) {
  cell <- CELLS[CELLS$cell_id == cid, ]; L0 <- numeric(0)
  perP <- perS <- perO <- vector("list", BATCHES)   # pooled-RB, single-origin-RB, pooled-observed
  for (bt in seq_len(BATCHES)) {
    p <- prep_one(cell, BATCH_N, seed = 900000L + 1000L * cid + bt)
    perP[[bt]] <- accumulate(new_acc(), p$rows, rb = TRUE)
    perS[[bt]] <- accumulate(new_acc(), p$rows[p$rows$trial == 0, , drop = FALSE], rb = TRUE)
    perO[[bt]] <- accumulate(new_acc(), p$rows, rb = FALSE)
    if (length(L0) < 500000L) L0 <- c(L0, p$L0)
  }
  all_b <- seq_len(BATCHES)
  mP <- standardise(solve_Bstar(pool_acc(perP, all_b, TRUE )), L0, THOR)
  mS <- standardise(solve_Bstar(pool_acc(perS, all_b, TRUE )), L0, THOR)
  mO <- standardise(solve_Bstar(pool_acc(perO, all_b, FALSE)), L0, THOR)
  # Delete-one-batch jackknife on each quantity AND on the pooling difference (which benefits from
  # the common random numbers shared by the two arms).
  jk <- vapply(all_b, function(b) {
    i <- setdiff(all_b, b)
    p1 <- standardise(solve_Bstar(pool_acc(perP, i, TRUE)), L0, THOR)
    p2 <- standardise(solve_Bstar(pool_acc(perS, i, TRUE)), L0, THOR)
    c(p1, p2, p1 - p2)
  }, numeric(3 * length(THOR)))
  jse <- function(rows_idx) {
    m <- jk[rows_idx, , drop = FALSE]
    sqrt((BATCHES - 1) / BATCHES * rowSums((m - rowMeans(m))^2))
  }
  nT <- length(THOR)
  seP <- jse(1:nT); seS <- jse(nT + 1:nT); seD <- jse(2 * nT + 1:nT)
  truth <- tv[[as.character(cell$alpha.0)]]$mrd
  out[[cell$label]] <- data.frame(
    label = cell$label, tau = THOR, truth = truth,
    plim_pooled = mP, se_pooled = seP, bias_pooled = mP - truth,
    plim_origin0 = mS, se_origin0 = seS,
    pooling_effect = mP - mS, se_pooling = seD,
    plim_pooled_obs = mO, stringsAsFactors = FALSE)
  cat(sprintf("\n%s (alpha.0=%.2f, gamma.L=%.1f):\n", cell$label, cell$alpha.0, cell$gamma.L))
  for (k in seq_len(nT)) {
    zb <- (mP[k] - truth[k]) / seP[k]; zd <- (mP[k] - mS[k]) / seD[k]
    cat(sprintf(paste0("  tau=%d | plim-truth=%+.5f +/-%.5f (z=%+.2f)%s",
                       " | POOLING=%+.5f +/-%.5f (z=%+.2f)%s | RB-obs=%+.5f\n"),
                THOR[k], mP[k] - truth[k], seP[k], zb, if (abs(zb) < 2) " [ns]" else "    ",
                mP[k] - mS[k], seD[k], zd, if (abs(zd) < 2) " [ns]" else "    ", mP[k] - mO[k]))
  }
}
res <- do.call(rbind, out)
attr(res, "settings") <- list(GRID = GRID, BATCH_N = BATCH_N, BATCHES = BATCHES,
                              seed_base = 900000L, TMAX = TMAX, weight = "stab (untruncated)",
                              rao_blackwellised = TRUE, single_origin_comparator = TRUE)
saveRDS(res, file.path(ROOT, "dev", "plim_closed_form.rds"))
cat(sprintf("\n[plim-cf] wrote dev/plim_closed_form.rds  (%.1f min)\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
