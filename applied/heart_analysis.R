# Applied example: Stanford heart transplant data (survival::jasa) analysed by
# sequential trial emulation with the steCI additive-hazards estimator + all six
# CI families. Treatment = heart transplant (time-varying, absorbing); outcome =
# death. The canonical immortal-time / time-varying-treatment dataset, fully
# public (ships in base R 'survival').
#
#   HEART_TIES=simultaneous Rscript applied/heart_analysis.R     (or HEART_TIES=timereg)
#
# HEART_TIES selects how tied event times are handled for EVERY printed estimate (see
# applied/heart_common.R); HEART_OUT overrides the output path (default
# applied/heart_results.rds). Backs Table heart, Table heart_mrd and the Section 5 prose.
#
# Reproducibility. Under HEART_TIES = "timereg" every timereg::aalen() fit on tied data
# consumes random numbers (timereg breaks ties by jitter), so the seed is set before the
# FIRST fit. (It used to sit after the point estimate and the printed analytic table, which
# were therefore irreproducible, and was itself dead code: ci_bootstrap(seed = 1) reseeds
# before its first draw.) ci_bootstrap() and ci_lef() still reseed with seed = 1; every other
# draw follows from heart_seed(). Under "simultaneous" only the resampling draws are random.
suppressMessages({library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang)})
for (f in list.files("steCI/R", pattern = "[.]R$", full.names = TRUE)) source(f)
source("applied/heart_common.R")
TIES <- heart_ties()
OUT <- Sys.getenv("HEART_OUT", unset = "applied/heart_results.rds")
heart_seed()                                         # before the first random draw
cat(sprintf("Tie handling: %s ; seed %d ; output %s\n", TIES, HEART_SEED, OUT))

K <- 6; W <- 30; TH <- 1:K
d <- heart_to_ste(w = W, K = K)
cat(sprintf("Person-period data: %d rows, %d patients, %d events; visits=%d (%d-day)\n",
            nrow(d), length(unique(d$id)), sum(d$event), K, W))
cat("Transplanted (ever) at baseline of each trial (A==1 share by visit):\n")
print(round(tapply(d$A, d$visit, mean), 3))

# ---- point estimates: computed ONCE, reported, and used to centre M4 and M5 ----
# time.grid = W: event times on the same trial-relative day are exact ties (not last-bit
# near-ties from time.stop - trial), in the point fit, every refit and the M6 fit alike.
est <- ste_estimate(d, n.visit = K, t.hor = TH, estimand = "both", ties = TIES, time.grid = W)
ds  <- attr(est, "dat.seq")
stopifnot(heart_near_ties(ds) == 0)
cat(sprintf("\nConstant hazard difference of transplant (per 30-day visit): %.4f\n", est$gamma))
cat("Marginal risk difference (death), always- vs never-transplant, by month:\n")
print(round(setNames(est$mrd, paste0("m", TH)), 4))

# ---- E1 hazard difference: M1 model / M2 robust / M3 cluster ----
analytic <- ci_analytic(ds, ties = TIES)
cat("\n== E1 constant hazard difference: analytic SEs ==\n")
print(analytic, digits = 3)

# ---- M4 bootstrap, M5 jackknife (n=103 feasible), M6 analytic+multiplier+band ----
bo <- ci_bootstrap(d, n.visit = K, t.hor = TH, B = 1000, estimand = "both", seed = 1,
                   ties = TIES, time.grid = W, point = est)
jk <- ci_jackknife(d, n.visit = K, t.hor = TH, estimand = "both", max_n = 200,
                   ties = TIES, time.grid = W, point = est)
lf <- ci_lef(d, n.visit = K, t.hor = TH, B = 1000, seed = 1, time.grid = W)

cat("\n== E1 hazard difference: bootstrap & jackknife ==\n")
print(bo$gamma, digits = 3); print(jk$gamma, digits = 3)
cat(sprintf("bootstrap: %d failed refits; %s refits skipped a singular risk set\n",
            bo$n_fail, format(bo$n_singular)))

cat("\n== E2 marginal risk difference (death) by month: bootstrap ==\n")
print(bo$mrd, digits = 3, row.names = FALSE)
cat("\n== E2 MRD: M6 analytic (row/cluster) + multiplier + band ==\n")
print(lf$mrd, digits = 3, row.names = FALSE)
print(lf$band, digits = 3, row.names = FALSE)

# ---- the two main-text tables, from the objects above ----
# Under "simultaneous" the MRD point estimate IS the M6 fit (same data, same simultaneous-tie
# fit), so the estimate column and the M6 columns come from one estimator; tolerance keyed to
# the magnitude of the estimates (last-bit differences in summation order only).
if (TIES == "simultaneous")
  stopifnot(max(abs(est$mrd - lf$point$mrd)) <= 1e3 * .Machine$double.eps * max(abs(est$mrd)))
se_of <- function(type) lf$mrd$se[lf$mrd$type == type]
tab_heart <- data.frame(M1_model = analytic$se[analytic$method == "model"],
                        M2_robust = analytic$se[analytic$method == "robust"],
                        M3_cluster = analytic$se[analytic$method == "cluster"],
                        M5_jackknife = jk$gamma$se, M4_bootstrap = bo$gamma$se[1])
tab_heart_mrd <- data.frame(horizon = TH, estimate = est$mrd,
                            se_row = se_of("analytic_row"), se_cluster = se_of("analytic_cluster"),
                            se_multiplier = se_of("multiplier"),
                            se_bootstrap = bo$mrd$se[bo$mrd$type == "percentile"])
cat(sprintf("\n== Table heart (caption: gamma-hat = %.3f per visit) ==\n", est$gamma))
print(round(tab_heart, 3), row.names = FALSE)
cat("\n== Table heart_mrd ==\n")
print(round(tab_heart_mrd, 3), row.names = FALSE)

saveRDS(list(d = d, est = est, analytic = analytic, boot = bo, jack = jk, lef = lf,
             K = K, W = W, ties = TIES, seed = HEART_SEED,
             tab_heart = tab_heart, tab_heart_mrd = tab_heart_mrd,
             versions = c(R = R.version.string,
                          timereg = as.character(utils::packageVersion("timereg")),
                          survival = as.character(utils::packageVersion("survival")))),
        OUT)
cat(sprintf("\nSaved %s\n", OUT))
