# Web Table S5 (appendix tab:stdvar) -- PERSIST PER-REPLICATE ESTIMATES.
#
# Re-runs the exact replicate stream of dev/experiment_stdvar.R (same 7 cells, same
# seeds 1000*ci + r, same dat_sim / ste_prepare / .aalen_if / multiplier Xi draws) but
# stores, per replicate x horizon, the ESTIMATES that coverage is computed from:
#   mrd, se_f (fixed / conditional), se_r (random / unconditional),
#   multiplier quantiles qf_lo, qf_hi (fixed) and qr_lo, qr_hi (random),
# plus the fitted cumulative coefficients. NO TRUTH IS READ HERE: the truth is an input
# to SCORING (score_s5.R), never to simulation.
#
# Package code is sourced from SNAP (env S5_SNAPSHOT; default steCI/R of this repository).
# The stored run used a frozen copy of the package code; the code paths it uses
# (dat_sim, ste_prepare, .aalen_if with
# its default singular = "error") compute the same numbers in the current package.
#
# Resumable: one .rds per (cell, rep-block); an existing chunk is skipped. Writes are
# atomic (tmp + rename). PSOCK cluster (mclapply does not fork on Windows).
#
#   Rscript dev/s5_rescore/run_s5_reps.R            # full run
#   S5_SMOKE=1 Rscript dev/s5_rescore/run_s5_reps.R # 2 reps per cell -> reps_smoke/
#   S5_WORKERS=6 ...                                # worker count (default 6)
suppressMessages(library(parallel))

ROOT  <- normalizePath(Sys.getenv("REPO_ROOT", "."))
SNAP  <- normalizePath(Sys.getenv("S5_SNAPSHOT", file.path(ROOT, "steCI/R")))
SMOKE <- identical(Sys.getenv("S5_SMOKE"), "1")
OUT   <- file.path(ROOT, "dev/s5_rescore", if (SMOKE) "reps_smoke" else "reps")
NW    <- as.integer(Sys.getenv("S5_WORKERS", "6"))
stopifnot(dir.exists(SNAP), length(list.files(SNAP, "[.]R$")) >= 10)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

TH <- 1:5
# representative grid -- IDENTICAL to dev/experiment_stdvar.R lines 23-28
cells <- data.frame(
  n       = c(300L, 300L, 300L, 1000L, 1000L, 1000L, 5000L),
  alpha.0 = c(0.15, 0.25, 0.40, 0.15,  0.25,  0.40,  0.25),
  reps    = c(400L, 400L, 400L, 400L,  400L,  400L,  200L),
  Bmult   = 200L, gamma.0 = -1, gamma.L = 1.5)
cells$ev <- c("low", "med", "high")[match(cells$alpha.0, c(0.15, 0.25, 0.40))]
if (SMOKE) cells$reps <- 2L

# task list: blocks of reps (smaller blocks at n=5000 so the load balances)
tasks <- list()
for (ci in seq_len(nrow(cells))) {
  bs <- if (cells$n[ci] >= 5000L) 5L else 20L
  st <- seq(1L, cells$reps[ci], by = bs)
  for (s in st) {
    e <- min(s + bs - 1L, cells$reps[ci])
    f <- file.path(OUT, sprintf("cell%d_r%03d-%03d.rds", ci, s, e))
    if (!file.exists(f)) tasks[[length(tasks) + 1L]] <- list(ci = ci, r0 = s, r1 = e, file = f)
  }
}
cat(sprintf("[%s] %d chunk(s) to run (%s), %d workers, snapshot %s\n",
            format(Sys.time()), length(tasks), if (SMOKE) "SMOKE" else "FULL", NW, Sys.getenv("S5_SNAPSHOT", "steCI/R")))
if (!length(tasks)) { cat("NOTHING TO DO\n"); quit(save = "no") }

# one replicate: the estimation half of experiment_stdvar.R::one(), verbatim arithmetic,
# with the coverage indicators replaced by the quantities they are computed from.
est_one <- function(d, Bmult, seed, TH) {
  ds <- ste_prepare(d); fit <- .aalen_if(ds, TH)
  L0 <- d$L[d$visit == 1]
  Lc <- ds$L.baseline[match(fit$clusters, as.character(ds$id))]
  K  <- dim(fit$IFc)[1]; nt <- length(TH); a <- 0.025
  set.seed(seed); Xi <- matrix(sample(c(-1, 1), K * Bmult, replace = TRUE), K, Bmult)
  out <- vector("list", nt)
  for (s in seq_len(nt)) {
    B0 <- fit$B[s, 1]; BA <- fit$B[s, 2]; BL <- fit$B[s, 3]
    S1 <- exp(-(B0 + BA + BL * L0)); S0 <- exp(-(B0 + BL * L0))
    mrd <- mean(S0 - S1)
    g   <- c(mean(S1 - S0), mean(S1), mean(L0 * (S1 - S0)))
    base <- as.vector(fit$IFc[, , s] %*% g)
    S1c <- exp(-(B0 + BA + BL * Lc)); S0c <- exp(-(B0 + BL * Lc))
    stdv <- ((S0c - S1c) - mean(S0c - S1c)) / K
    rand <- base + stdv
    se_f <- sqrt(sum(base^2)); se_r <- sqrt(sum(rand^2))
    mf <- mrd + as.vector(crossprod(base, Xi))
    mr <- mrd + as.vector(crossprod(rand, Xi))
    qf <- stats::quantile(mf, c(a, 1 - a), names = FALSE)
    qr <- stats::quantile(mr, c(a, 1 - a), names = FALSE)
    out[[s]] <- data.frame(tau = TH[s], B0 = B0, BA = BA, BL = BL, mrd = mrd,
                           se_f = se_f, se_r = se_r, qf_lo = qf[1], qf_hi = qf[2],
                           qr_lo = qr[1], qr_hi = qr[2], K = K, n_L0 = length(L0),
                           n_rows = nrow(ds))
  }
  do.call(rbind, out)
}

run_task <- function(tk, cells, TH) {
  cell <- cells[tk$ci, ]
  rows <- list()
  for (r in tk$r0:tk$r1) {
    seed <- 1000L * tk$ci + r
    key <- data.frame(ci = tk$ci, n = cell$n, alpha0 = cell$alpha.0, ev = cell$ev,
                      rep = r, seed = seed)
    status <- "ok"; msg <- NA_character_; est <- NULL
    d <- tryCatch(dat_sim(n = cell$n, alpha.0 = cell$alpha.0, gamma.0 = cell$gamma.0,
                          gamma.L = cell$gamma.L, seed = seed),
                  error = function(e) { msg <<- conditionMessage(e); NULL })
    if (is.null(d)) status <- "dat_sim_error" else {
      est <- tryCatch(est_one(d, cell$Bmult, seed, TH),
                      error = function(e) { msg <<- conditionMessage(e); NULL })
      if (is.null(est)) status <- "est_error"
    }
    if (is.null(est))
      est <- data.frame(tau = TH, B0 = NA_real_, BA = NA_real_, BL = NA_real_, mrd = NA_real_,
                        se_f = NA_real_, se_r = NA_real_, qf_lo = NA_real_, qf_hi = NA_real_,
                        qr_lo = NA_real_, qr_hi = NA_real_, K = NA_integer_,
                        n_L0 = NA_integer_, n_rows = NA_integer_)
    rows[[length(rows) + 1L]] <- cbind(key, status = status, msg = msg, est,
                                       stringsAsFactors = FALSE)
  }
  res <- do.call(rbind, rows)
  attr(res, "meta") <- list(R = R.version.string, timereg = as.character(packageVersion("timereg")),
                            time = format(Sys.time()), host = Sys.info()[["nodename"]])
  tmp <- paste0(tk$file, ".tmp")
  saveRDS(res, tmp); file.rename(tmp, tk$file)
  sprintf("%s ok=%d/%d", basename(tk$file), sum(res$status == "ok") / length(TH),
          tk$r1 - tk$r0 + 1L)
}

cl <- makePSOCKcluster(NW)
on.exit(stopCluster(cl), add = TRUE)
clusterExport(cl, c("SNAP", "est_one", "run_task"))
invisible(clusterEvalQ(cl, {
  suppressMessages({library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang)})
  for (f in list.files(SNAP, pattern = "[.]R$", full.names = TRUE)) source(f)
  NULL
}))
t0 <- Sys.time()
# load-balanced; log each chunk as it lands (via clusterApplyLB returning in order is
# fine -- the per-chunk files are the durable record)
res <- clusterApplyLB(cl, tasks, function(tk, cells, TH) {
  tryCatch(run_task(tk, cells, TH), error = function(e) paste("CHUNK ERROR", tk$file, conditionMessage(e)))
}, cells = cells, TH = TH)
cat(unlist(res), sep = "\n")
cat(sprintf("[%s] elapsed %.1f min\n", format(Sys.time()),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
n_done <- length(list.files(OUT, "^cell.*[.]rds$"))
cat(sprintf("chunk files present: %d\n", n_done))
cat("RUN_S5_REPS DONE\n")
