# Aggregate all intermediate per-task RDS files into replicate-level and
# metric-level summaries.
#   Rscript simulation/aggregate.R

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
getc <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) x else d }
REPO_ROOT <- getc("REPO_ROOT", ".")

suppressMessages({ library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang) })
load_steCI <- function(root) {
  if (requireNamespace("steCI", quietly = TRUE)) { library(steCI); return("installed") }
  cands <- unique(c(file.path(root, "steCI", "R"),
                    file.path(root, "..", "steCI", "R"),
                    "steCI/R", "../steCI/R"))
  hit <- cands[dir.exists(cands)][1]
  if (is.na(hit)) stop(
    "Cannot find the steCI package. Place the 'steCI/' folder at REPO_ROOT (",
    normalizePath(root, mustWork = FALSE), ") or as a sibling of 'simulation/'.",
    call. = FALSE)
  for (f in list.files(hit, pattern = "[.]R$", full.names = TRUE)) source(f)
  hit
}
load_steCI(REPO_ROOT)
source(file.path(REPO_ROOT, "simulation", "scenarios.R"))
source(file.path(REPO_ROOT, "simulation", "metrics.R"))

OUT <- getc("SIM_OUT", file.path(REPO_ROOT, "simulation", "output", "intermediate"))
RES <- getc("SIM_RESULTS", file.path(REPO_ROOT, "simulation", "output", "results"))
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

files <- list.files(OUT, pattern = "^sim_.*\\.rds$", full.names = TRUE)
if (!length(files)) stop("no intermediate files in ", OUT)
cat("[aggregate] reading", length(files), "files\n")

reps <- do.call(rbind, lapply(files, function(f) readRDS(f)$reps))
saveRDS(reps, file.path(RES, "replicates_all.rds"))
cat("[aggregate] total replicate rows:", nrow(reps), "\n")

# Score against the exact closed-form truths used by the paper. Historical
# diagnostic helpers can select a different existing cache with SIM_TRUTH_CACHE.
# A missing cache is an input error; do not silently substitute a Monte Carlo draw.
cache_file <- getc("SIM_TRUTH_CACHE", file.path(REPO_ROOT, "dev", "truth_cache_exact.rds"))
if (!file.exists(cache_file))
  stop("truth cache not found: ", cache_file,
       "\nRebuild the exact cache with Rscript dev/truth_closed_form.R, or set SIM_TRUTH_CACHE.")
tv <- readRDS(cache_file)
cat("[aggregate] using cached truths:", cache_file, "\n")
truth_fun <- truth_lookup(tv)
m <- summarize_metrics(reps, truth_fun)

saveRDS(m, file.path(RES, "metrics.rds"))
utils::write.csv(m$ci, file.path(RES, "metrics_ci.csv"), row.names = FALSE)
utils::write.csv(m$empirical, file.path(RES, "metrics_empirical.csv"), row.names = FALSE)
# Simultaneous (across-tau) coverage -- the band's claim is a simultaneous one, so it must be
# measured, not inferred from pointwise averages. NULL when there is only one horizon.
if (!is.null(m$joint)) utils::write.csv(m$joint, file.path(RES, "metrics_joint.csv"), row.names = FALSE)
cat("[aggregate] wrote metrics to", RES, "\n")

# Replicate-level outputs that failed definitional validity are guarded out of the SCALE
# statistics (emp_sd/mean_se/se_ratio) but never silently dropped -- surface the count here so a
# re-run cannot hide it. See the note above .valid_estimate() in metrics.R.
n_bad_est <- sum(m$empirical$n_invalid_est); n_bad_se <- sum(m$ci$n_invalid_se)
if (n_bad_est || n_bad_se) {
  cat("[aggregate] WARNING: guarded out", n_bad_est, "point estimates outside the parameter space and",
      n_bad_se, "invalid SEs (per method-cell). Raw values retained in the *_raw columns.\n")
}

# Headline: coverage of the hazard-difference CIs, averaged over scenarios.
hd <- m$ci[m$ci$estimand == "hd", ]
agg <- aggregate(cbind(coverage, se_ratio, mean_width) ~ method, data = hd, FUN = mean)
cat("\n== Hazard-difference CI coverage (mean over scenarios) ==\n")
print(agg)
