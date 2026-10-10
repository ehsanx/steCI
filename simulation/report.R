# Reproduce Tables 3-5 and the coverage summaries from the shipped exact-truth
# metrics, without rerunning the simulation or requiring raw replicate files.
#   Rscript simulation/report.R
#
# Env: REPO_ROOT (default .), RESULTS_DIR (default data/exact/results).
# Set SIM_OUT explicitly to score per-task raw chunks instead. That path uses
# dev/truth_cache_exact.rds unless SIM_TRUTH_CACHE names another cache explicitly;
# a missing cache is an error, never a request for new Monte Carlo truths.
ROOT <- Sys.getenv("REPO_ROOT", ".")
source(file.path(ROOT, "simulation/scenarios.R"))
scen <- make_scenarios()

OUT <- Sys.getenv("SIM_OUT", "")
if (nzchar(OUT)) {
  source(file.path(ROOT, "simulation/metrics.R"))
  files <- list.files(OUT, pattern = "^sim_.*[.]rds$", full.names = TRUE)
  if (!length(files)) stop("no intermediate files in ", OUT)
  reps <- do.call(rbind, lapply(files, function(f) readRDS(f)$reps))
  if (is.null(reps) || !nrow(reps)) stop("no replicate rows in ", OUT)
  truth_cache <- Sys.getenv("SIM_TRUTH_CACHE",
                           file.path(ROOT, "dev/truth_cache_exact.rds"))
  if (!file.exists(truth_cache))
    stop("truth cache not found: ", truth_cache,
         "\nRebuild the exact cache with Rscript dev/truth_closed_form.R.")
  tv <- readRDS(truth_cache)
  needed <- as.character(sort(unique(reps$alpha.0)))
  if (!all(needed %in% names(tv))) stop("truth cache lacks an input alpha.0 level")
  cat("[report] scoring raw chunks against:", truth_cache, "\n")
  m <- summarize_metrics(reps, truth_lookup(tv))
  cat(sprintf("Loaded %d files, %d rows, %d scenarios (n: %s)\n",
              length(files), nrow(reps), length(unique(reps$scenario)),
              paste(sort(unique(reps$n)), collapse = ",")))
  cat("Raw input can be a partial run; check its completed replicate counts.\n")
} else {
  RES <- Sys.getenv("RESULTS_DIR", file.path(ROOT, "data/exact/results"))
  metrics_file <- file.path(RES, "metrics.rds")
  if (!file.exists(metrics_file)) stop("metrics file not found: ", metrics_file)
  m <- readRDS(metrics_file)
  cat("[report] reading shipped metrics:", metrics_file, "\n")
}
if (!is.list(m) || !all(c("ci", "empirical") %in% names(m)))
  stop("metrics must contain ci and empirical data frames")

add_design <- function(x) {
  needed <- setdiff(c("n", "event_rate"), names(x))
  if (!length(needed)) return(x)
  if (!all(x$scenario %in% scen$scenario)) stop("metrics contain unknown scenarios")
  merge(x, scen[, c("scenario", needed), drop = FALSE], by = "scenario")
}
ci <- add_design(m$ci)
emp <- add_design(m$empirical)
e1_methods <- c("model", "robust", "cluster", "boot_percentile")
e2_methods <- c("m6_analytic_row", "m6_analytic_cluster", "m6_multiplier",
                "m6_band", "m6a2_analytic_cluster", "m6a2_multiplier",
                "boot_percentile")
hd <- ci[ci$estimand == "hd" & ci$method %in% e1_methods, ]
mrd <- ci[ci$estimand == "mrd" & ci$method %in% e2_methods, ]

cat("\n== E1 marginal constant hazard difference: bias & empirical SD ==\n")
print(aggregate(cbind(bias, emp_sd) ~ n, emp[emp$estimand == "hd", ], mean),
      digits = 3, row.names = FALSE)

cat("\n== Table 3: E1 coverage / bias-eliminated coverage / SE ratio by method x n ==\n")
table3 <- aggregate(cbind(coverage, coverage_be, se_ratio) ~ method + n, hd, mean)
print(table3, digits = 3, row.names = FALSE)

cat("\n== Table 4: E1 coverage by method x event rate x n ==\n")
table4 <- aggregate(coverage ~ method + event_rate + n, hd, mean)
print(table4, digits = 3, row.names = FALSE)

cat("\n== Table 5: E2 MRD coverage by method x n, averaged over horizons ==\n")
table5 <- aggregate(coverage ~ method + n, mrd, mean)
print(table5, digits = 3, row.names = FALSE)

cat("\n== E2 MRD coverage by method x horizon x n ==\n")
print(aggregate(coverage ~ method + tau + n, mrd, mean),
      digits = 3, row.names = FALSE)

cat("\n== E2 MRD coverage by method x horizon, averaged over scenarios ==\n")
print(aggregate(coverage ~ method + tau, mrd, mean),
      digits = 3, row.names = FALSE)

if (!is.null(m$joint) && nrow(m$joint)) {
  joint <- m$joint[m$joint$estimand == "mrd", ]
  joint$hits <- joint$joint_coverage * joint$n_curves
  simultaneous <- aggregate(cbind(hits, n_curves) ~ method, joint, sum)
  simultaneous$coverage <- simultaneous$hits / simultaneous$n_curves
  cat("\n== E2 simultaneous coverage, weighted by complete curves ==\n")
  print(simultaneous[, c("method", "coverage", "n_curves")],
        digits = 4, row.names = FALSE)
}
