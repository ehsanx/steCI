# M8 apply step. Builds a NEW truth cache in which gamma_true is the multi-seed mean from
# dev/truth_multiseed.rds, leaving every E2 survival truth byte-identical, then re-aggregates
# into a SIBLING results directory so nothing canonical is overwritten until the diff is read.
#
#   Rscript dev/apply_m8_truth.R            # build cache + aggregate into results_m8/
#
# Promotion to canonical is a deliberate, separate step -- see dev/compare_m8.R first.
P      <- "data"
CANON  <- file.path(P, "truth_cache.rds")
# This historical comparison uses its own cache explicitly, rather than the
# published exact-truth cache used by aggregate.R by default.
NEWDIR <- file.path(P, "m8")
NEWRES <- file.path(NEWDIR, "results")

stopifnot(file.exists(CANON), file.exists("dev/truth_multiseed.rds"))
ms <- readRDS("dev/truth_multiseed.rds")
tc <- readRDS(CANON)

dir.create(NEWDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(NEWRES, showWarnings = FALSE, recursive = TRUE)

new <- tc
for (a0 in names(tc)) {
  stopifnot(!is.null(ms[[a0]]))
  # E2 truths (t.hor, surv1, surv0, survdiff, mrd) are carried over UNTOUCHED.
  new[[a0]]$gamma_true      <- ms[[a0]]$mean
  new[[a0]]$gamma_mcse      <- ms[[a0]]$mcse   # the SE M8 asks to be reported
  new[[a0]]$gamma_sd_1draw  <- ms[[a0]]$sd
  new[[a0]]$gamma_nseed     <- length(ms[[a0]]$gamma)
  new[[a0]]$gamma_cached_old<- tc[[a0]]$gamma_true
}

# Assert the E2 truths really are unchanged -- this is the whole basis for M7/M4 being
# independent of M8, so it is checked rather than trusted.
for (a0 in names(tc)) for (f in c("t.hor","surv1","surv0","survdiff","mrd"))
  stopifnot(identical(tc[[a0]][[f]], new[[a0]][[f]]))
cat("[m8] E2 survival truths verified identical across all", length(tc), "alpha.0 levels\n")

saveRDS(new, file.path(NEWDIR, "truth_cache.rds"))
cat("[m8] wrote", file.path(NEWDIR, "truth_cache.rds"), "\n")
for (a0 in names(new)) cat(sprintf("     a0=%-5s gamma %.6f -> %.6f  (shift %+.6f, MCSE %.6f)\n",
  a0, new[[a0]]$gamma_cached_old, new[[a0]]$gamma_true,
  new[[a0]]$gamma_true - new[[a0]]$gamma_cached_old, new[[a0]]$gamma_mcse))

cat("\n[m8] re-aggregating into", NEWRES, "(no canonical file is touched)\n")
stopifnot(file.exists(file.path(dirname(NEWRES), "truth_cache.rds")))  # aggregate must find OURS
Sys.setenv(REPO_ROOT = ".",
           SIM_OUT = file.path(P, "intermediate"),
           SIM_RESULTS = NEWRES,
           SIM_TRUTH_CACHE = file.path(NEWDIR, "truth_cache.rds"))
source("simulation/aggregate.R")
cat("\n[m8] done. Canonical results untouched. Next: Rscript dev/compare_m8.R\n")
