# What actually moved when the single-frozen-draw truth defect is fixed. Compares the canonical
# metrics against the multi-seed-truth re-aggregation, and answers the question:
# does the E1 "residual is largely point-estimator bias" story survive?
#
#   Rscript dev/compare_m8.R
P   <- "data"
OLD <- file.path(P, "results")
NEW <- file.path(P, "m8", "results")
stopifnot(dir.exists(OLD), dir.exists(NEW))
nof <- function(s) as.integer(sub("^n([0-9]+)_.*$", "\\1", s))

rd <- function(d, f) read.csv(file.path(d, f), stringsAsFactors = FALSE)
oc <- rd(OLD, "metrics_ci.csv");        nc <- rd(NEW, "metrics_ci.csv")
oe <- rd(OLD, "metrics_empirical.csv"); ne <- rd(NEW, "metrics_empirical.csv")
for (d in list(oc, nc)) stopifnot(nrow(d) > 0)
stopifnot(nrow(oc) == nrow(nc), nrow(oe) == nrow(ne))

cat("=========== 1. E2 MUST BE UNCHANGED (the independence claim) ===========\n")
m_o <- oc[oc$estimand == "mrd", ]; m_n <- nc[nc$estimand == "mrd", ]
cat(sprintf("  max |coverage diff| over %d E2 cells: %.3e\n", nrow(m_o),
            max(abs(m_o$coverage - m_n$coverage))))
cat(sprintf("  max |truth diff|    over %d E2 cells: %.3e\n", nrow(m_o),
            max(abs(m_o$truth - m_n$truth))))
cat("  -> if either is non-zero, the E2 results are NOT independent of the truth fix.\n")

cat("\n=========== 2. E1 truth ===========\n")
h_o <- oc[oc$estimand == "hd", ]; h_n <- nc[nc$estimand == "hd", ]
cat(sprintf("  old truths: %s\n", paste(sprintf("%.6f", sort(unique(h_o$truth))), collapse = ", ")))
cat(sprintf("  new truths: %s\n", paste(sprintf("%.6f", sort(unique(h_n$truth))), collapse = ", ")))

cat("\n=========== 3. E1 bias, by sample size ===========\n")
eo <- oe[oe$estimand == "hd", ]; en <- ne[ne$estimand == "hd", ]
eo$nn <- nof(eo$scenario); en$nn <- nof(en$scenario)
cmp <- merge(aggregate(bias ~ nn, eo, mean), aggregate(bias ~ nn, en, mean),
             by = "nn", suffixes = c("_old", "_new"))
cmp$pct_change <- 100 * (abs(cmp$bias_new) - abs(cmp$bias_old)) / abs(cmp$bias_old)
print(cmp, digits = 4, row.names = FALSE)

cat("\n=========== 4. THE DECIDING NUMBER: does the E1 bias story survive? ===========\n")
cat("recentring gain = coverage_be - coverage, averaged over scenarios, by method x n.\n")
cat("The paper leans on E1 gains of 0.016-0.054 being LARGE (bias explains E1's residual),\n")
cat("against E2 <=0.006 (variance, unexplained). If the E1 gains collapse, that story fails.\n\n")
gain <- function(d) { d$nn <- nof(d$scenario); d$g <- d$coverage_be - d$coverage
  aggregate(g ~ method + nn, d[d$estimand == "hd", ], mean) }
go <- gain(oc); gn <- gain(nc)
g <- merge(go, gn, by = c("method", "nn"), suffixes = c("_old", "_new"))
g <- g[order(g$method, g$nn), ]
print(g, digits = 3, row.names = FALSE)
cat(sprintf("\n  E1 recentring gain, RANGE over all method x n:  old %.4f-%.4f   new %.4f-%.4f\n",
            min(go$g), max(go$g), min(gn$g), max(gn$g)))
e2g <- function(d) { d$g <- d$coverage_be - d$coverage; max(aggregate(g ~ method, d[d$estimand=="mrd",], mean)$g) }
cat(sprintf("  E2 recentring gain, max over methods:          %.4f (unchanged by construction)\n", e2g(nc)))

cat("\n=========== 5. E1 raw coverage (truth-dependent) ===========\n")
cv <- merge(aggregate(coverage ~ method + nn, transform(h_o, nn = nof(scenario)), mean),
            aggregate(coverage ~ method + nn, transform(h_n, nn = nof(scenario)), mean),
            by = c("method", "nn"), suffixes = c("_old", "_new"))
cv$diff <- cv$coverage_new - cv$coverage_old
print(cv[order(cv$method, cv$nn), ], digits = 3, row.names = FALSE)

cat("\n=========== 6. TRUTH-FREE quantities MUST NOT MOVE ===========\n")
for (col in c("se_ratio", "coverage_be")) {
  a <- oc[[col]]; b <- nc[[col]]
  cat(sprintf("  %-12s max |diff| = %.3e %s\n", col, max(abs(a - b), na.rm = TRUE),
              if (max(abs(a - b), na.rm = TRUE) < 1e-12) "OK" else "<-- UNEXPECTED"))
}
cat("\n(se_ratio and coverage_be use emp_sd/mean_est, never truth. Any movement here is a bug.)\n")
