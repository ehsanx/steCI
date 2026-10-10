# BETWEEN-STREAM uncertainty of the closed-form plim at long horizons (worst cell).
#
# WHY. dev/plim_closed_form.rds reports a delete-one-batch jackknife SE (se_pooled ~ 0.0016 at
# tau=5). That is a WITHIN-stream measure and it UNDERSTATES the true uncertainty: two independent
# seed streams already gave tau=5 = +0.00700 +/- 0.00157 (seed 900000) and +0.00121 +/- 0.00458
# (seed 700000) -- their jackknife SEs disagree 3x and the point estimates are 3.7x the smaller SE
# apart. The heavy-tailed IPACW weights mean the long-horizon estimate is dominated by rare
# heavy-weight batches, so whether a stream HITS one is a between-stream effect the jackknife cannot
# see. This script runs several INDEPENDENT streams and reports the between-stream SD directly, so
# the honest tau>=3 uncertainty is a measured number rather than an anecdote.
#
#   Rscript dev/plim_seedcheck.R
# Env: PLIM_SC_STREAMS (5), PLIM_SC_BATCHES (25), PLIM_SC_BATCH_N (20000), PLIM_SC_GRID (4000).
suppressMessages({library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang)})
ROOT <- "."
for (f in list.files(file.path(ROOT, "steCI", "R"), pattern = "[.]R$", full.names = TRUE)) source(f)
source(file.path(ROOT, "simulation", "scenarios.R"))
tv <- readRDS(file.path(ROOT, "dev", "truth_cache_exact.rds"))

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
TMAX <- 5; GRID <- geti("PLIM_SC_GRID", 4000L); DT <- TMAX / GRID; THOR <- T_HOR
stopifnot(GRID %% TMAX == 0L)
A0 <- 0.15; GL <- 3.0; G0 <- 0.0; AA <- DGM_FIXED$alpha.A; AL <- DGM_FIXED$alpha.L; AU <- DGM_FIXED$alpha.U
TRUTH <- tv[["0.15"]]$mrd

scadd <- function(len, idx, val) {
  keep <- idx >= 1L & idx <= len & is.finite(val); if (!any(keep)) return(numeric(len))
  r <- rowsum(val[keep], idx[keep]); out <- numeric(len); out[as.integer(rownames(r))] <- r[, 1]; out }
new_acc <- function() list(gram = replicate(6, numeric(GRID + 1L), simplify = FALSE),
                           dlam = replicate(3, numeric(GRID + 1L), simplify = FALSE))
accumulate <- function(acc, rows) {
  a <- rows$a; b <- rows$b; w <- rows$w; X <- cbind(1, rows$A, rows$L)
  jlo <- pmax(1L, floor(a/DT)+1L); jhi <- pmin(GRID, floor(b/DT)); ok <- jhi >= jlo
  ent <- list(c(1,1),c(1,2),c(1,3),c(2,2),c(2,3),c(3,3))
  for (k in seq_along(ent)) { p<-ent[[k]][1]; q<-ent[[k]][2]; v<-w*X[,p]*X[,q]
    acc$gram[[k]] <- acc$gram[[k]] + scadd(GRID+1L, jlo[ok], v[ok]) - scadd(GRID+1L, jhi[ok]+1L, v[ok]) }
  for (p in 1:3) { v <- w*X[,p]*rows$h*DT
    acc$dlam[[p]] <- acc$dlam[[p]] + scadd(GRID+1L, jlo[ok], v[ok]) - scadd(GRID+1L, jhi[ok]+1L, v[ok]) }
  acc }
solve_Bstar <- function(acc) {
  gc <- lapply(acc$gram, function(x) cumsum(x)[1:GRID]); dc <- lapply(acc$dlam, function(x) cumsum(x)[1:GRID])
  dB <- matrix(0, GRID, 3)
  for (j in 1:GRID) { G <- matrix(c(gc[[1]][j],gc[[2]][j],gc[[3]][j],gc[[2]][j],gc[[4]][j],gc[[5]][j],gc[[3]][j],gc[[5]][j],gc[[6]][j]),3,3)
    dB[j,] <- tryCatch(solve(G, c(dc[[1]][j],dc[[2]][j],dc[[3]][j])), error=function(e) c(NA,NA,NA)) }
  dB[!is.finite(dB)] <- 0; apply(dB, 2, cumsum) }
standardise <- function(B, L0, thor) { Jt <- pmin(GRID, floor(thor/DT))
  vapply(seq_along(thor), function(k) { b0<-B[Jt[k],1]; bA<-B[Jt[k],2]; bL<-B[Jt[k],3]
    mean(exp(-(b0+bL*L0))) - mean(exp(-(b0+bA+bL*L0))) }, numeric(1)) }
pool_acc <- function(per, idx) { a <- new_acc()
  for (k in 1:6) a$gram[[k]] <- Reduce(`+`, lapply(per[idx], function(z) z$gram[[k]]))
  for (p in 1:3) a$dlam[[p]] <- Reduce(`+`, lapply(per[idx], function(z) z$dlam[[p]])); a }
prep <- function(n, seed) {
  d <- dat_sim(n=n, n.visit=DGM_FIXED$n.visit, gamma.0=G0, gamma.L=GL, alpha.0=A0, alpha.A=AA,
               alpha.L=AL, alpha.U=AU, seed=seed)
  set.seed(seed); U <- stats::rnorm(n, 0, 0.1); ds <- ste_prepare(d, n.visit=DGM_FIXED$n.visit)
  hpos <- pmax(A0 + AA*ds$A + AL*ds$L + AU*U[ds$id], 0)
  list(rows=data.frame(a=ds$time.new, b=ds$time.stop.new, A=ds$A.baseline, L=ds$L.baseline, w=ds$ipw.stab.cum, h=hpos),
       L0 = d$L[d$visit==1]) }

NSTR <- geti("PLIM_SC_STREAMS", 5L); B_B <- geti("PLIM_SC_BATCHES", 25L); B_N <- geti("PLIM_SC_BATCH_N", 20000L)
# Distinct seed bases, none equal to the 900000/700000 already run (so these ADD independent streams).
BASES <- c(110000L, 220000L, 330000L, 440000L, 550000L, 660000L, 880000L)[seq_len(NSTR)]
cat(sprintf("=== between-stream check, cell 1 | %d streams x %d batches x n=%d | GRID=%d ===\n", NSTR, B_B, B_N, GRID))
M <- matrix(NA_real_, NSTR, length(THOR)); JK <- matrix(NA_real_, NSTR, length(THOR)); t0 <- Sys.time()
for (s in seq_len(NSTR)) {
  base <- BASES[s]; L0 <- numeric(0); per <- vector("list", B_B)
  for (bt in seq_len(B_B)) { p <- prep(B_N, seed = base + bt); per[[bt]] <- accumulate(new_acc(), p$rows)
    if (length(L0) < 400000L) L0 <- c(L0, p$L0) }
  ab <- seq_len(B_B); m <- standardise(solve_Bstar(pool_acc(per, ab)), L0, THOR) - TRUTH
  jk <- vapply(ab, function(b) standardise(solve_Bstar(pool_acc(per, setdiff(ab,b))), L0, THOR) - TRUTH, numeric(length(THOR)))
  se <- sqrt((B_B-1)/B_B * rowSums((jk - rowMeans(jk))^2))
  M[s,] <- m; JK[s,] <- se
  cat(sprintf("  stream %d (base %d): tau5 bias=%+.5f (jack %.5f)\n", s, base, m[5], se[5]))
}
# Add the two existing independent streams for the pooled summary.
M2 <- rbind(M, c(0.00700, NA, NA, NA, NA)); M2[nrow(M2), ] <- c(-1.09e-5, 9.71e-5, 9.33e-4, 3.380e-3, 7.003e-3)
M2 <- rbind(M2, c(-9.20e-6, 4.82e-4, 1.735e-3, 2.612e-3, 1.208e-3))   # seed 700000 (ablation baseline)
cat("\n=== SUMMARY over", NSTR, "new +2 existing =", NSTR+2, "independent streams ===\n")
for (k in seq_along(THOR)) {
  vals <- M2[, k]; bsd <- sd(vals); bmean <- mean(vals); bse <- bsd/sqrt(length(vals))
  meanjk <- mean(c(JK[, k]), na.rm = TRUE)
  cat(sprintf("  tau=%d  mean %+.5f | BETWEEN-stream SD %.5f (SE of mean %.5f) | mean within-stream jackknife %.5f | understatement %.1fx\n",
              THOR[k], bmean, bsd, bse, meanjk, bsd/meanjk))
}
saveRDS(list(new=M, jackknife=JK, bases=BASES, all=M2, thor=THOR,
             settings=list(GRID=GRID, BATCH_N=B_N, BATCHES=B_B, streams=NSTR)),
        file.path(ROOT, "dev", "plim_seedcheck.rds"))
cat(sprintf("\n[seedcheck] wrote dev/plim_seedcheck.rds (%.1f min)\n", as.numeric(difftime(Sys.time(), t0, units="mins"))))
