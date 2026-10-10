# Does pinning max.clust = NULL change the PUBLISHED M3 numbers at n=5000?
# The canonical run was produced WITHOUT the pin, so if the pin moves M3 there, the code no
# longer reproduces the paper -- which would be a worse defect than the one being fixed.
suppressMessages({library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang)})
for (f in list.files("steCI/R", pattern = "[.]R$", full.names = TRUE)) source(f)
source("simulation/scenarios.R")

SEED_BASE <- 20250228L; N_REPS <- 1000L
set.seed(SEED_BASE); all_seeds <- sample.int(.Machine$integer.max, N_REPS)
g <- make_scenarios(); cell <- g[g$scenario == "n5000_evlow_cfweak_pvmed", ]

set.seed(all_seeds[1])
d <- dat_sim(n = cell$n, n.visit = DGM_FIXED$n.visit, gamma.0 = cell$gamma.0,
             gamma.L = cell$gamma.L, alpha.0 = cell$alpha.0, alpha.A = DGM_FIXED$alpha.A,
             alpha.L = DGM_FIXED$alpha.L, alpha.U = DGM_FIXED$alpha.U)
est <- ste_estimate(d, n.visit = DGM_FIXED$n.visit, t.hor = T_HOR,
                    weight = "stab", estimand = "both")
ds <- attr(est, "dat.seq")
w  <- ds$ipw.stab.cum
cat(sprintf("n=%d, stacked rows %d, subjects %d\n", cell$n, nrow(ds), length(unique(ds$id))))

fit <- function(mc) {
  a <- if (is.null(mc))
    timereg::aalen(survival::Surv(time.new, time.stop.new, event) ~
                     const(A.baseline) + L.baseline, data = ds, n.sim = 0,
                   weights = w, robust = 1, clusters = ds$id, max.clust = NULL)
  else
    timereg::aalen(survival::Surv(time.new, time.stop.new, event) ~
                     const(A.baseline) + L.baseline, data = ds, n.sim = 0,
                   weights = w, robust = 1, clusters = ds$id)
  sqrt(a$robvar.gamma[1, 1])
}
pinned <- fit(NULL); default <- fit("default")
cat(sprintf("\n  M3 SE with max.clust = NULL (pinned) : %.10f\n", pinned))
cat(sprintf("  M3 SE with timereg default          : %.10f\n", default))
cat(sprintf("  difference                          : %.3e  (rel %.3e)\n",
            pinned - default, abs(pinned - default) / default))

can <- readRDS("data/exact/results/replicates_all.rds")
cs <- can$se[can$scenario == "n5000_evlow_cfweak_pvmed" & can$method == "cluster" &
             can$estimand == "hd" & can$rep == 1]
cat(sprintf("\n  canonical published value           : %.17g\n", cs))
cat(sprintf("  pinned                              : %.17g\n", pinned))
cat(sprintf("  default                             : %.17g\n", default))

# This block previously used isTRUE(all.equal(default, cs, tolerance = 1e-10)) and printed
# "default matches published? FALSE" -- a FALSE ALARM, and not a numerical one: all.equal
# also compares ATTRIBUTES, and sqrt(a$robvar.gamma[1, 1]) carries a dimname that cs does
# not, so all.equal returned an attribute-mismatch string and isTRUE() read it as failure
# while the values agreed to 6.5 ulp. Compare numerics numerically.
#
# Key the tolerance to MAGNITUDE, never to exact zero.
# Identical arithmetic in a different evaluation order (BLAS threading, chunk order)
# differs in the last few ulp; the floor is |x|*eps, not 0. We allow 1e4 ulp, which is
# ~2.2e-12 relative here -- far tighter than the 1.6e-2 effect being tested, so it
# still fails loudly if the two fits genuinely disagree.
ulp <- function(x) abs(x) * .Machine$double.eps
rel <- function(a, b) abs(a - b) / abs(b)
cat(sprintf("  |default - published| = %.3e  (%.1f ulp, rel %.3e)\n",
            abs(default - cs), abs(default - cs) / ulp(cs), rel(default, cs)))
cat(sprintf("  |pinned  - published| = %.3e  (%.1f ulp, rel %.3e)\n",
            abs(pinned - cs), abs(pinned - cs) / ulp(cs), rel(pinned, cs)))
cat(sprintf("\n  default reproduces published (<=1e4 ulp)? %s\n",
            abs(default - cs) <= 1e4 * ulp(cs)))
cat(sprintf("  pinned differs from published, as expected? %s\n",
            abs(pinned - cs) > 1e4 * ulp(cs)))
