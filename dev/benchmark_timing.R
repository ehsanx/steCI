# Per-operation timing that backs the computational-cost claims in the paper
# (M6 costs one influence-function fit; M4 refits the whole pipeline B times).
#   Rscript dev/benchmark_timing.R [n]        (default n = 5000)
args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args)) as.integer(args[1]) else 5000L
suppressMessages({library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang)})
for (f in list.files("steCI/R", pattern = "[.]R$", full.names = TRUE)) source(f)

secs <- c()
tm <- function(lbl, e) {
  t <- Sys.time(); v <- force(e)
  s <- as.numeric(difftime(Sys.time(), t, units = "secs")); secs[[lbl]] <<- s
  cat(sprintf("  %-32s %7.2f s\n", lbl, s))
  invisible(v)
}

cat(sprintf("Timing at n = %d\n", n))
d  <- dat_sim(n = n, seed = 99)
ds <- tm("ste_prepare (once)", ste_prepare(d))
one_prep_fit <- tm("ste_estimate both (prep + 2 fits)", ste_estimate(d, estimand = "both"))
tm("aalen_if (M6 fit incl. IFs)", .aalen_if(ds, 1:5))
tm("ci_lef (M6), B = 200", ci_lef(d, B = 200, seed = 1))
refit <- secs[["ste_estimate both (prep + 2 fits)"]]; m6 <- secs[["ci_lef (M6), B = 200"]]
cat(sprintf(paste0("\nM4 (nonparametric bootstrap) per replicate ~ B x (prep + 2 fits) = 200 x %.1f s = %.0f min (extrapolated);\n",
                   "M6 (ci_lef, B = 200) per replicate = %.1f s (measured). Stored runs: dev/benchmark_timing.log.\n"),
            refit, 200 * refit / 60, m6))
