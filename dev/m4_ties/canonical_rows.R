# Full-precision values of every PUBLISHED row the M4 tie decision rule refers to, recomputed from
# the canonical exact aggregation (data/exact/results). Output: canonical_rows.rds (+ printout).
# The printed manuscript values are reproduced here so the rounding-boundary trigger (c) can be
# evaluated against full precision, not against the rounded print.
RES <- Sys.getenv("M4_CANON_RESULTS",
  "data/exact/results")
stopifnot(dir.exists(RES), basename(dirname(RES)) == "exact")
m <- readRDS(file.path(RES, "metrics.rds"))
ci <- m$ci; jt <- m$joint
sc <- do.call(rbind, strsplit(ci$scenario, "_"))
ci$n <- as.integer(sub("n", "", sc[, 1])); ci$ev <- sub("ev", "", sc[, 2])
out <- list()
# Table e1 (hd), rows M1-M4: coverage, coverage_be, se_ratio by n (mean over 27 scenarios)
hd <- ci[ci$estimand == "hd", ]
out$tab_e1 <- aggregate(cbind(coverage, coverage_be, se_ratio) ~ method + n, data = hd, FUN = mean)
# Table e1_eventrate: coverage by n x event rate (9 scenarios)
out$tab_e1_eventrate <- aggregate(coverage ~ method + n + ev, data = hd, FUN = mean)
# Table e2 (mrd): coverage by n (mean over 27 scenarios x 5 horizons)
mr <- ci[ci$estimand == "mrd", ]
out$tab_e2 <- aggregate(coverage ~ method + n, data = mr, FUN = mean)
# Web Table S3: coverage by tau over 81 scenarios
out$tab_S3 <- aggregate(coverage ~ method + tau, data = mr, FUN = mean)
# Web Table S8: se_ratio by tau over 81 scenarios
out$tab_S8 <- aggregate(se_ratio ~ method + tau, data = mr, FUN = mean)
# simultaneous coverage, weighted by n_curves, overall and by n
jt$n <- as.integer(sub("n", "", sub("_.*", "", jt$scenario)))
wmean <- function(d) sum(d$joint_coverage * d$n_curves) / sum(d$n_curves)
out$simult <- do.call(rbind, lapply(split(jt, jt$method), function(d)
  data.frame(method = d$method[1], all = wmean(d), n300 = wmean(d[d$n == 300, ]),
             n1000 = wmean(d[d$n == 1000, ]), n5000 = wmean(d[d$n == 5000, ]),
             n_curves = sum(d$n_curves))))
out$simult <- out$simult[order(-out$simult$all), ]
saveRDS(out, file.path(Sys.getenv("M4_OUT", "dev/m4_ties"), "canonical_rows.rds"))
op <- options(width = 160, digits = 6)
b <- function(d) d[grepl("boot|cluster|multiplier|m6", d$method), ]
cat("== Table e1 (hd) ==\n"); print(out$tab_e1)
cat("== Table e1_eventrate (hd), M3/M4 ==\n"); print(out$tab_e1_eventrate[out$tab_e1_eventrate$method %in% c("cluster", "boot_percentile", "boot_normal", "boot_basic"), ])
cat("== Table e2 (mrd) ==\n"); print(out$tab_e2)
cat("== Web Table S3 ==\n"); print(reshape(out$tab_S3, idvar = "method", timevar = "tau", direction = "wide"))
cat("== Web Table S8 ==\n"); print(reshape(out$tab_S8, idvar = "method", timevar = "tau", direction = "wide"))
cat("== Simultaneous (weighted by n_curves) ==\n"); print(out$simult)
options(op)
