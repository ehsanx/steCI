# Re-aggregate the canonical results against the EXACT closed-form truths (dev/truth_closed_form.R).
#
# Writes to a SIBLING results directory. Nothing canonical is touched until the diff has been read
# -- this changes every E1 number and some E2 ones, so it is not a change to make in place.
#
# This historical raw-chunk promotion helper selects the new cache explicitly.
# Rebuilding the shipped exact truths alone requires dev/truth_closed_form.R.
#
#   Rscript dev/apply_exact_truth.R
P      <- "data"
NEWDIR <- file.path(P, "exact")
NEWRES <- file.path(NEWDIR, "results")
stopifnot(file.exists("dev/truth_cache_exact.rds"))
dir.create(NEWRES, showWarnings = FALSE, recursive = TRUE)

new <- readRDS("dev/truth_cache_exact.rds")
old <- readRDS(file.path(P, "truth_cache.rds"))

# The exact cache must be shape-compatible with what truth_lookup() expects, or metrics silently
# score against the wrong thing. Check the contract rather than assume it.
for (k in names(old)) {
  stopifnot(!is.null(new[[k]]), identical(new[[k]]$t.hor, old[[k]]$t.hor),
            length(new[[k]]$mrd) == length(old[[k]]$mrd),
            is.finite(new[[k]]$gamma_true))
}
cat("[exact] cache contract OK for", length(new), "alpha.0 levels\n")
for (k in names(new)) cat(sprintf("   a0=%-5s gamma %+.8f (was %+.8f, shift %+.6f)\n",
                                  k, new[[k]]$gamma_true, old[[k]]$gamma_true,
                                  new[[k]]$gamma_true - old[[k]]$gamma_true))

saveRDS(new, file.path(NEWDIR, "truth_cache.rds"))
stopifnot(file.exists(file.path(dirname(NEWRES), "truth_cache.rds")))

cat("\n[exact] re-aggregating into", NEWRES, "-- canonical results untouched\n")
Sys.setenv(REPO_ROOT = ".", SIM_OUT = file.path(P, "intermediate"),
           SIM_RESULTS = NEWRES,
           SIM_TRUTH_CACHE = file.path(NEWDIR, "truth_cache.rds"))
source("simulation/aggregate.R")
cat("\n[exact] done. Next: Rscript dev/compare_exact.R\n")
