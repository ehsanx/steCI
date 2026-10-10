# dev/heart_grid_sensitivity.R
# Backs the visit-grid sensitivity of the applied heart example:
#   - Web Table S9 (appendix, tab:heartgrid): the 6-rung analytic ladder
#     (w = 30/20/15/12/10/7 days; 6-month MRD, row/cluster SEs, ESS, max weight).
#   - Section 5 main-text sentences: MRD attenuates -0.09 -> ~+0.02 toward the null;
#     row-level SE 26% below cluster at 30-day, 45% below at weekly; the finest-rung
#     jackknife/bootstrap SEs corroborate the analytic clustered SE.
#
# Message: the effect estimate is grid-dependent (a within-interval treatment-timing
# coarsening effect, NOT weight degeneracy -- weights stay healthy, ESS >96%), while
# the row-vs-cluster SE understatement the example illustrates is invariant to the
# grid and, if anything, widens. gamma is a per-visit quantity and is NOT
# magnitude-comparable across w, so the ladder is anchored on the horizon-fixed MRD.
#
#   HEART_TIES=simultaneous Rscript dev/heart_grid_sensitivity.R     (or HEART_TIES=timereg)
#
# HEART_TIES is THE tie switch shared with applied/heart_analysis.R (applied/heart_common.R):
# every estimate below -- the ladder's MRD and SEs and the corroborating refits -- uses it, so
# S9 and Tables heart / heart_mrd come from one estimator. HEART_GRID_OUT overrides the output
# path (default dev/heart_grid_sensitivity_results.rds).
# Deterministic: heart_seed() precedes the first fit (under "timereg" every fit on tied data
# draws jitter), and the parallel bootstrap/jackknife are reproducible for a fixed worker
# count via clusterSetRNGStream (NW = 8 PSOCK workers). The w = 30 rung repeats
# applied/heart_analysis.R's first two fits in the same order after the same seed, so it
# reproduces that script's MRD and M2/M3 SEs exactly under either setting.
if (!file.exists("applied/heart_common.R")) {        # not launched from the repo root
  f <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])
  setwd(dirname(dirname(normalizePath(f))))
}
suppressMessages({library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang); library(parallel)})
for (f in list.files("steCI/R", pattern = "[.]R$", full.names = TRUE)) source(f)
source("applied/heart_common.R")
TIES <- heart_ties()
OUT <- Sys.getenv("HEART_GRID_OUT", unset = "dev/heart_grid_sensitivity_results.rds")
NW <- 8L
heart_seed()                                         # before the first random draw
cat(sprintf("Tie handling: %s ; seed %d ; %d workers ; output %s\n", TIES, HEART_SEED, NW, OUT))

## ---- 1. Analytic ladder (Web Table S9) ----
fit_rung <- function(w) {
  K <- ceiling(180 / w); TH <- 1:K; d <- suppressMessages(heart_to_ste(w, K))
  est <- ste_estimate(d, n.visit = K, t.hor = TH, estimand = "both", ties = TIES, time.grid = w)
  ds <- attr(est, "dat.seq"); stopifnot(heart_near_ties(ds) == 0)
  ana <- ci_analytic(ds, ties = TIES)
  list(w = w, K = K, TH = TH, d = d, gamma = est$gamma, mrd = est$mrd, ana = ana,
       wc = ds$ipw.stab.cum)
}
rungs <- lapply(c(30, 20, 15, 12, 10, 7), fit_rung)
ladder <- do.call(rbind, lapply(rungs, function(r) {
  j6 <- which.min(abs(r$TH * r$w - 180)); se <- r$ana$se; wc <- r$wc
  data.frame(w = r$w, K = r$K, mrd6 = r$mrd[j6], se_row = se[2], se_cluster = se[3],
             row_below_clu_pct = 100 * (1 - se[2] / se[3]),
             ess_pct = (sum(wc)^2 / sum(wc^2)) / length(wc) * 100, max_w = max(wc),
             gamma = r$gamma)
}))
cat("=== Web Table S9: analytic ladder ===\n")
print(transform(ladder, mrd6 = round(mrd6, 3), se_row = round(se_row, 3),
                se_cluster = round(se_cluster, 3), row_below_clu_pct = round(row_below_clu_pct, 1),
                ess_pct = round(ess_pct, 1), max_w = round(max_w, 2), gamma = NULL),
      row.names = FALSE)
cat("--- as printed in Web Table S9 ---\n")
print(data.frame(w = ladder$w, K = ladder$K, mrd6 = sprintf("%+.2f", ladder$mrd6),
                 se_row = sprintf("%.3f", ladder$se_row), se_cluster = sprintf("%.3f", ladder$se_cluster),
                 row_below = sprintf("%.0f%%", ladder$row_below_clu_pct),
                 ess = sprintf("%.0f%%", ladder$ess_pct), max_wt = sprintf("%.1f", ladder$max_w)),
      row.names = FALSE)

## ---- 2. Finest-rung jackknife + bootstrap corroboration ----
root <- getwd()
cl <- makeCluster(NW)
clusterExport(cl, c("root", "TIES"))
invisible(clusterEvalQ(cl, { setwd(root)
  suppressMessages({library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang)})
  for (f in list.files("steCI/R", pattern = "[.]R$", full.names = TRUE)) source(f); NULL }))
clusterSetRNGStream(cl, HEART_SEED)
corrob <- function(r) {                              # reuses the ladder's fit of this rung
  d <- r$d; K <- r$K; TH <- r$TH; w <- r$w
  clusterExport(cl, c("d", "K", "TH", "w"), envir = environment())
  B <- 600; ch <- rep(B %/% NW, NW); ch[seq_len(B %% NW)] <- ch[seq_len(B %% NW)] + 1
  gb <- unlist(parLapply(cl, ch, function(nb) { rs <- .make_resampler(d); g <- rep(NA, nb)
    for (b in seq_len(nb)) { e <- tryCatch(ste_estimate(rs(), n.visit = K, t.hor = TH, estimand = "hd",
                                                        ties = TIES, time.grid = w),
                                           error = function(e) NULL)
      if (!is.null(e)) g[b] <- e$gamma }; g }))
  ids <- unique(d$id)
  gj <- unlist(parLapply(cl, ids, function(k) { e <- tryCatch(ste_estimate(d[d$id != k, ], n.visit = K,
      t.hor = TH, estimand = "hd", ties = TIES, time.grid = w), error = function(e) NULL)
    if (is.null(e)) NA else e$gamma }))
  nj <- sum(!is.na(gj))
  list(summary = data.frame(w = w, se_row = r$ana$se[2], se_cluster = r$ana$se[3],
         se_jackknife = sqrt((nj - 1) / nj * sum((gj - mean(gj, na.rm = TRUE))^2, na.rm = TRUE)),
         se_bootstrap = sd(gb, na.rm = TRUE),
         # M5 as ci_jackknife() defines it (centred on the full-sample estimate), for reference
         se_jackknife_pointcentred = sqrt((nj - 1) / nj * sum((gj - r$gamma)^2, na.rm = TRUE)),
         n_boot_fail = sum(is.na(gb)), n_jack_fail = sum(is.na(gj))),
       boot = gb, jack = gj)
}
cor <- lapply(rungs[match(c(30, 10, 7), sapply(rungs, `[[`, "w"))], corrob)
stopCluster(cl)
corrob_tab <- do.call(rbind, lapply(cor, `[[`, "summary"))
cat("\n=== finest-rung corroboration (row vs subject-level resampling SEs) ===\n")
print(transform(corrob_tab, se_row = round(se_row, 3), se_cluster = round(se_cluster, 3),
                se_jackknife = round(se_jackknife, 3), se_bootstrap = round(se_bootstrap, 3),
                se_jackknife_pointcentred = round(se_jackknife_pointcentred, 3)), row.names = FALSE)

saveRDS(list(ties = TIES, seed = HEART_SEED, NW = NW, ladder = ladder, corrob = corrob_tab,
             boot_draws = setNames(lapply(cor, `[[`, "boot"), c("w30", "w10", "w7")),
             jack_values = setNames(lapply(cor, `[[`, "jack"), c("w30", "w10", "w7")),
             rungs = lapply(rungs, function(r) r[c("w", "K", "gamma", "mrd", "ana")]),
             versions = c(R = R.version.string,
                          timereg = as.character(utils::packageVersion("timereg")),
                          survival = as.character(utils::packageVersion("survival")))),
        OUT)
cat(sprintf("\nSaved %s\n", OUT))
