# Extract, ONCE, the canonical rows of the main simulation the M4 tie experiment pairs against:
#   canon_boot.rds  : every point / boot_* row (all 81 scenarios x all reps), jitter arm of record
#   canon_nfail.rds : per (scen_id, rep) canonical boot_n_fail (used for the advance-mode check)
# Source: the canonical exact aggregation (data/exact/results; replicates_all.rds is not in git). The boot_* and
# point columns are truth-free, so the exact vs results/ distinction does not bite here, but the
# canonical path is used regardless.
RES <- Sys.getenv("M4_CANON_RESULTS",
  "data/exact/results")
stopifnot(dir.exists(RES), basename(dirname(RES)) == "exact",
          file.exists(file.path(RES, "replicates_all.rds")))
r <- readRDS(file.path(RES, "replicates_all.rds"))
keep <- r$method %in% c("point", "boot_percentile", "boot_basic", "boot_normal")
cb <- r[keep, c("scen_id", "scenario", "rep", "n", "alpha.0", "gamma.0", "gamma.L", "estimand",
                "tau", "method", "estimate", "se", "lower", "upper", "boot_n_fail", "failed_methods")]
rownames(cb) <- NULL
OUT <- Sys.getenv("M4_DATA", "dev/m4_ties")
saveRDS(cb, file.path(OUT, "canon_boot.rds"))
nf <- cb[cb$method == "boot_percentile" & cb$estimand == "hd", c("scen_id", "rep", "boot_n_fail")]
saveRDS(nf, file.path(OUT, "canon_nfail.rds"))
cat("canon_boot rows:", nrow(cb), " canon_nfail rows:", nrow(nf), "\n")
