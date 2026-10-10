# Calibrate the additive-hazard DGM: report the negative-hazard truncation
# fraction and the achieved event rate at each (event-rate x prevalence x
# confounding) combination, at the extremes that stress hazard non-negativity.
#
#   Rscript simulation/calibrate_dgm.R
#
# Use this to confirm event-rate levels before launching the full simulation.
# A neg_haz fraction above ~0.02-0.03 means the additive hazard is being
# clipped hard enough to distort the estimand -- raise alpha.0 for that level.

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
getc <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) x else d }
REPO_ROOT <- getc("REPO_ROOT", ".")
suppressMessages({ library(dplyr); library(tidyr); library(rlang) })
src <- unique(c(file.path(REPO_ROOT, "steCI", "R"), file.path(REPO_ROOT, "..", "steCI", "R"),
                "steCI/R", "../steCI/R"))
src <- src[dir.exists(src)][1]
for (f in list.files(src, pattern = "[.]R$", full.names = TRUE)) source(f)

n <- geti("CAL_N", 4000L)
grid <- expand.grid(alpha.0 = c(0.15, 0.25, 0.40),
                    gamma.0 = c(0.0, -1.0, -2.0),   # high/med/low prevalence
                    gamma.L = c(0.5, 1.5, 3.0))     # weak/mod/strong confounding

cat(sprintf("DGM calibration (n=%d per cell)\n", n))
cat(sprintf("%-8s %-8s %-8s %-10s %-10s %-10s\n",
            "alpha0", "gamma0", "gammaL", "neg_haz", "event_rt", "prev_A1"))
for (i in seq_len(nrow(grid))) {
  d <- dat_sim(n = n, gamma.0 = grid$gamma.0[i], gamma.L = grid$gamma.L[i],
               alpha.0 = grid$alpha.0[i], seed = 1000 + i)
  nid <- length(unique(d$id))
  er  <- sum(d$event) / nid
  prev <- mean(d$A[d$visit == 1])   # baseline treatment prevalence
  cat(sprintf("%-8.2f %-8.0f %-8.1f %-10.3f %-10.3f %-10.3f\n",
              grid$alpha.0[i], grid$gamma.0[i], grid$gamma.L[i],
              attr(d, "neg_haz_frac"), er, prev))
}
