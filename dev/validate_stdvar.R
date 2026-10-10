# Validate the random-standardisation (unconditional) MRD variance option in
# ci_lef(): the analytic cluster and multiplier SEs under standardization="random"
# should EXCEED the fixed (conditional) ones and closely match the nonparametric
# bootstrap SE, which resamples subjects and so already captures the sampling
# variability of the empirical standardisation population.
#   Rscript dev/validate_stdvar.R
suppressMessages({library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang)})
for (f in list.files("steCI/R", pattern = "[.]R$", full.names = TRUE)) source(f)

set.seed(11)
d  <- dat_sim(n = 800, alpha.0 = 0.25, gamma.0 = -1, gamma.L = 1.5, seed = 11)
fx <- ci_lef(d, B = 500, standardization = "fixed",  seed = 1)
rd <- ci_lef(d, B = 500, standardization = "random", seed = 1)
bt <- ci_bootstrap(d, B = 300, estimand = "mrd", t.hor = 1:5, seed = 7)

getse <- function(o, type) o$mrd$se[o$mrd$type == type]
cat("horizon                :", sprintf("%7d", 1:5), "\n")
cat("cluster  SE (fixed)    :", sprintf("%.4f", getse(fx, "analytic_cluster")), "\n")
cat("cluster  SE (random)   :", sprintf("%.4f", getse(rd, "analytic_cluster")), "\n")
cat("multiplier SE (fixed)  :", sprintf("%.4f", getse(fx, "multiplier")), "\n")
cat("multiplier SE (random) :", sprintf("%.4f", getse(rd, "multiplier")), "\n")
cat("bootstrap SE (M4 ref)  :", sprintf("%.4f", bt$mrd$se[bt$mrd$type == "normal"]), "\n")
cat("\nratio random/bootstrap :", sprintf("%.3f", getse(rd,"analytic_cluster")/bt$mrd$se[bt$mrd$type=="normal"]), "\n")
cat("ratio fixed /bootstrap :", sprintf("%.3f", getse(fx,"analytic_cluster")/bt$mrd$se[bt$mrd$type=="normal"]), "\n")
cat("\nVALIDATE_STDVAR DONE\n")
