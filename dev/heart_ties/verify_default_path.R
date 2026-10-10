# Does the tie switch / time.grid change anything for the continuous-time SIMULATION?
# It must not: the package defaults (ties = "timereg", time.grid = NULL) are the published path.
#
# Reruns complete canonical replicates -- EVERY method the main simulation recorded (point,
# M1-M3, M4 percentile/basic/normal at B = 200, M6 A1 row/cluster/multiplier/band, M6 A2) --
# through simulation/one_replicate.R with the steCI/R code in STECI_R_DIR (default: this
# repo's steCI/R), and compares row by row with the canonical per-replicate record
#   data/exact/results/replicates_all.rds
# (91 MB, not in git; see the README, "Large files"; CANON_REPLICATES overrides). Seeds as in run_sim.R: set.seed(20250228); sample.int(.Machine$integer.max, 1000).
#
#   Rscript dev/heart_ties/verify_default_path.R <scenario> <rep> [<scenario> <rep> ...]
#   STECI_R_DIR=<other steCI/R copy> Rscript dev/heart_ties/verify_default_path.R ...
#
# Tolerance keyed to magnitude: a row reproduces if every
# |difference| <= 1e-12 -- far above the last-ulp floor of
# BLAS/platform reordering (~1e-16 * |x|), far below any printed digit.
suppressMessages({library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang)})
RDIR <- Sys.getenv("STECI_R_DIR", unset = "steCI/R")
for (f in list.files(RDIR, pattern = "[.]R$", full.names = TRUE)) source(f)
source("simulation/scenarios.R"); source("simulation/one_replicate.R")
CANON <- Sys.getenv("CANON_REPLICATES", "data/exact/results/replicates_all.rds")
if (!file.exists(CANON)) stop("canonical replicate record not found: ", CANON)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2L, length(args) %% 2L == 0L)
jobs <- data.frame(scenario = args[c(TRUE, FALSE)], rep = as.integer(args[c(FALSE, TRUE)]))

set.seed(20250228L); all_seeds <- sample.int(.Machine$integer.max, 1000L)
canon <- readRDS(CANON)
g <- make_scenarios()
key <- function(x) paste(x$estimand, ifelse(is.na(x$tau), "NA", sprintf("%.6f", x$tau)), x$method)
TOL <- 1e-12
cat(sprintf("steCI code: %s\n", RDIR))
all_ok <- TRUE
for (i in seq_len(nrow(jobs))) {
  sc <- jobs$scenario[i]; rp <- jobs$rep[i]
  cell <- g[g$scenario == sc, ]; stopifnot(nrow(cell) == 1L)
  cv <- canon[canon$scenario == sc & canon$rep == rp, ]
  t0 <- Sys.time()
  new <- run_one(cell, rep_seed = all_seeds[rp], B = 200L, do_jackknife = FALSE)
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  kc <- key(cv); kn <- key(new)
  stopifnot(!anyDuplicated(kc), !anyDuplicated(kn))
  common <- intersect(kc, kn)
  m <- match(common, kc); mn <- match(common, kn)
  cols <- c("estimate", "se", "lower", "upper")
  dif <- sapply(cols, function(cc) {
    a <- cv[[cc]][m]; b <- new[[cc]][mn]
    both_na <- is.na(a) & is.na(b)
    if (any(xor(is.na(a), is.na(b)))) return(Inf)
    max(c(0, abs(a - b)[!both_na]))
  })
  ok <- length(common) == nrow(cv) && length(common) == nrow(new) && all(dif <= TOL) &&
        identical(cv$failed_methods[1], new$failed_methods[1]) &&
        identical(as.integer(cv$boot_n_fail[1]), as.integer(new$boot_n_fail[1]))
  all_ok <- all_ok && ok
  cat(sprintf("%-32s rep %4d: %3d/%3d rows matched; max|diff| est %.1e se %.1e lower %.1e upper %.1e; failed_methods '%s'/'%s'; %.1f min -> %s\n",
              sc, rp, length(common), nrow(cv), dif["estimate"], dif["se"], dif["lower"], dif["upper"],
              cv$failed_methods[1], new$failed_methods[1], el, if (ok) "REPRODUCES" else "DOES NOT REPRODUCE"))
}
cat(if (all_ok) "\nALL REPLICATES REPRODUCE (<= 1e-12)\n" else "\nMISMATCH -- investigate before trusting the change\n")
