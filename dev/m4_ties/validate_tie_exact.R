# Validation of the M4 tie-experiment machinery (run BEFORE the pilot and the production run).
#
#   V0  ste_prepare_fast() is IDENTICAL (identical(), every column the fits read) to the snapshot's
#       ste_prepare(), on original data and on bootstrap resamples, n = 300 / 1000 / 5000.
#   V1  On TIE-FREE data (the original simulated datasets), the tie-exact fit (E1) equals timereg:
#       gamma, intZHZ, intZHdN, B(tau) and MRD(tau), across all 27 n=300 cells, 9 n=1000 cells and
#       3 n=5000 cells (seed = canonical replicate 1000's seed). Tolerance keyed to magnitude.
#   V2  On bootstrap resamples (ties present), E1 (from-scratch simultaneous) equals E2 (timereg on
#       the multiplicity-collapsed, tie-free data), and E2 has no residual ties.
#   V3  Timereg's jitter is live on the same resamples: the jittered fit differs from E1/E2.
#
#   Rscript dev/m4_ties/validate_tie_exact.R   ->  dev/m4_ties/validation.rds + printed summary
suppressMessages({ library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang) })
ROOT <- Sys.getenv("REPO_ROOT", ".")
HERE <- file.path(ROOT, "dev/m4_ties")
SNAP <- Sys.getenv("M4_SNAPSHOT", ROOT)    # code tree for steCI/R + simulation/scenarios.R
for (f in list.files(file.path(SNAP, "steCI/R"), "[.]R$", full.names = TRUE)) source(f)
source(file.path(SNAP, "simulation/scenarios.R"))
source(file.path(HERE, "tie_exact.R")); source(file.path(HERE, "prepare_fast.R"))

scen <- make_scenarios()
set.seed(20250228L); all_seeds <- sample.int(.Machine$integer.max, 1000L)
seed_v <- all_seeds[1000]                         # a replicate NOT used by pilot/production
cells <- rbind(scen[scen$n == 300L, ],
               scen[scen$n == 1000L & scen$gamma.0 == -1, ],
               scen[scen$n == 5000L & scen$gamma.0 == -1 & scen$gamma.L == 1.5, ])
TH <- 1:5
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300))
res <- list()
for (i in seq_len(nrow(cells))) {
  cl <- cells[i, ]
  set.seed(seed_v)
  d <- dat_sim(n = cl$n, n.visit = 5, gamma.0 = cl$gamma.0, gamma.L = cl$gamma.L,
               alpha.0 = cl$alpha.0, alpha.A = DGM_FIXED$alpha.A, alpha.L = DGM_FIXED$alpha.L,
               alpha.U = DGM_FIXED$alpha.U)
  ds0 <- ste_prepare(d); dsf <- ste_prepare_fast(d)
  v0_orig <- prep_identical(dsf, ds0)
  ti <- tie_info(ds0)
  # V1: tie-free original data, timereg vs E1
  tr <- ste_fit(ds0, d, t.hor = TH, estimand = "both")
  e1 <- fit_exact_simul(ds0, d, t.hor = TH)
  cum <- tr$fit_mrd$cum
  Btr <- t(vapply(TH, function(h) { k <- findInterval(h, cum[, 1]); cum[k, 2:4] }, numeric(3)))
  v1 <- c(gamma_abs = abs(e1$gamma - tr$gamma),
          gamma_rel = rel(e1$gamma, tr$gamma),
          intZHZ_rel = rel(e1$intZHZ, tr$fit_hd$intZHZ[1, 1]),
          intZHdN_rel = rel(e1$intZHdN, tr$fit_hd$intZHdN[1, 1]),
          B_absmax = max(abs(e1$B - Btr)), B_relmax = rel(e1$B, Btr),
          mrd_absmax = max(abs(e1$mrd - tr$mrd)))
  # V2/V3: two bootstrap resamples of this dataset
  rs <- make_resampler_idx(d)
  v2 <- NULL
  for (b in 1:2) {
    dl <- rs()
    dsb <- ste_prepare(dl); dsbf <- ste_prepare_fast(dl)
    e1b <- fit_exact_simul(dsb, dl, t.hor = TH)
    e2b <- fit_exact_collapse(dsb, dl, attr(dl, "orig_id"), t.hor = TH)
    trb <- ste_fit(dsb, dl, t.hor = TH, estimand = "both")
    tib <- tie_info(dsb)
    v2 <- rbind(v2, c(v0_boot = prep_identical(dsbf, dsb),
                      nties = unname(tib["nties"]), dt = unname(tib["dt"]),
                      ties_left_E2 = e2b$ties_left,
                      E1E2_gamma_abs = abs(e1b$gamma - e2b$gamma),
                      E1E2_intZHZ_rel = rel(e1b$intZHZ, e2b$intZHZ),
                      E1E2_mrd_absmax = max(abs(e1b$mrd - e2b$mrd)),
                      jit_gamma_abs = abs(trb$gamma - e2b$gamma),
                      jit_mrd_absmax = max(abs(trb$mrd - e2b$mrd))))
  }
  res[[i]] <- list(cell = cl$scenario, n = cl$n, v0_orig = v0_orig, orig_nties = unname(ti["nties"]),
                   v1 = v1, v2 = v2)
  cat(sprintf("%-28s V0 %s/%s  ties0=%d  V1 gamma_rel %.1e  B_rel %.1e  mrd_abs %.1e | V2 E1-E2 mrd %.1e gam %.1e | jitter mrd %.1e\n",
              cl$scenario, v0_orig, all(v2[, "v0_boot"] == 1), ti["nties"], v1["gamma_rel"],
              v1["B_relmax"], v1["mrd_absmax"], max(v2[, "E1E2_mrd_absmax"]),
              max(v2[, "E1E2_gamma_abs"]), max(v2[, "jit_mrd_absmax"])))
}
saveRDS(res, file.path(Sys.getenv("M4_OUT", HERE), "validation.rds"))
V1 <- do.call(rbind, lapply(res, `[[`, "v1")); V2 <- do.call(rbind, lapply(res, `[[`, "v2"))
cat("\n==== SUMMARY ====\n")
cat(sprintf("V0 ste_prepare_fast identical: original %d/%d, resamples %d/%d\n",
            sum(sapply(res, `[[`, "v0_orig")), length(res), sum(V2[, "v0_boot"]), nrow(V2)))
cat(sprintf("original datasets with event ties: %d of %d\n", sum(sapply(res, `[[`, "orig_nties") > 0), length(res)))
cat("V1 (tie-free, E1 vs timereg) maxima:\n"); print(signif(apply(V1, 2, max), 3))
cat("V2 (resamples, E1 vs E2) maxima and residual ties:\n")
print(signif(apply(V2[, c("E1E2_gamma_abs", "E1E2_intZHZ_rel", "E1E2_mrd_absmax", "ties_left_E2")], 2, max), 3))
cat("V3 jitter live (timereg-jitter minus E2): median |gamma| ", signif(median(V2[, "jit_gamma_abs"]), 3),
    " median max|mrd| ", signif(median(V2[, "jit_mrd_absmax"]), 3), " ; nties range ",
    paste(range(V2[, "nties"]), collapse = "-"), "\n")
