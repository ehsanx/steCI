# What moved when the Monte-Carlo truths were replaced by exact ones.
#   Rscript dev/compare_exact.R
P   <- "data"
OLD <- file.path(P, "results"); NEW <- file.path(P, "exact", "results")
stopifnot(dir.exists(OLD), dir.exists(NEW))
nof <- function(s) as.integer(sub("^n([0-9]+)_.*$", "\\1", s))
rd  <- function(d, f) read.csv(file.path(d, f), stringsAsFactors = FALSE)

oc <- rd(OLD,"metrics_ci.csv");        nc <- rd(NEW,"metrics_ci.csv")
oe <- rd(OLD,"metrics_empirical.csv"); ne <- rd(NEW,"metrics_empirical.csv")
stopifnot(nrow(oc)==nrow(nc), nrow(oe)==nrow(ne))

cat("=========== 1. TRUTH-FREE quantities MUST NOT MOVE ===========\n")
cat("se_ratio and coverage_be use emp_sd/mean_est, never truth. Movement here is a BUG.\n")
for (col in c("se_ratio","coverage_be","se_ratio_fold")) {
  if (!col %in% names(oc)) next
  d <- max(abs(oc[[col]] - nc[[col]]), na.rm = TRUE)
  cat(sprintf("  %-14s max|diff| = %.3e  %s\n", col, d, if (d < 1e-12) "OK" else "<-- UNEXPECTED"))
}

cat("\n=========== 2. E1: bias and raw coverage ===========\n")
eo <- oe[oe$estimand=="hd",]; en <- ne[ne$estimand=="hd",]
eo$nn <- nof(eo$scenario); en$nn <- nof(en$scenario)
b <- merge(aggregate(bias~nn,eo,mean), aggregate(bias~nn,en,mean), by="nn", suffixes=c("_MC","_exact"))
b$pct <- 100*(abs(b$bias_exact)-abs(b$bias_MC))/abs(b$bias_MC)
cat("mean E1 bias by n:\n"); print(b, digits=4, row.names=FALSE)

ho <- oc[oc$estimand=="hd",]; hn <- nc[nc$estimand=="hd",]
ho$nn <- nof(ho$scenario); hn$nn <- nof(hn$scenario)
cv <- merge(aggregate(coverage~method+nn,ho,mean), aggregate(coverage~method+nn,hn,mean),
            by=c("method","nn"), suffixes=c("_MC","_exact"))
cv$diff <- cv$coverage_exact - cv$coverage_MC
cat("\nE1 coverage by method x n:\n"); print(cv[order(cv$method,cv$nn),], digits=3, row.names=FALSE)

cat("\n=========== 3. THE DECIDING NUMBER for M5 ===========\n")
cat("Recentring gain = coverage_be - coverage. The paper argues the E1 residual is LARGELY BIAS\n")
cat("on the strength of E1 gains of 0.016-0.054 against E2's <=0.008. If the E1 gains collapse,\n")
cat("that argument weakens and M5 must be rewritten around the new numbers.\n\n")
gain <- function(d, es) { d$nn <- nof(d$scenario); d$g <- d$coverage_be - d$coverage
  aggregate(g~method+nn, d[d$estimand==es,], mean) }
for (es in c("hd","mrd")) {
  go <- gain(oc,es); gn <- gain(nc,es)
  cat(sprintf("  %s: range over method x n  --  MC truth %.4f to %.4f | EXACT truth %.4f to %.4f\n",
              toupper(es), min(go$g), max(go$g), min(gn$g), max(gn$g)))
}
g1o <- gain(oc,"hd"); g1n <- gain(nc,"hd")
g <- merge(g1o, g1n, by=c("method","nn"), suffixes=c("_MC","_exact"))
cat("\nE1 recentring gain, per method x n:\n"); print(g[order(g$method,g$nn),], digits=3, row.names=FALSE)

cat("\n=========== 4. E2: does anything move? ===========\n")
mo <- oe[oe$estimand=="mrd",]; mn <- ne[ne$estimand=="mrd",]
cat(sprintf("  max |E2 bias diff| over %d cells: %.3e\n", nrow(mo), max(abs(mo$bias-mn$bias))))
m5o <- oc[oc$estimand=="mrd",]; m5n <- nc[nc$estimand=="mrd",]
cat(sprintf("  max |E2 coverage diff|:            %.3e\n", max(abs(m5o$coverage-m5n$coverage))))
cat(sprintf("  mean |E2 coverage diff|:           %.3e\n", mean(abs(m5o$coverage-m5n$coverage))))

cat("\n=========== 5. G90: the four stale table entries ===========\n")
cat("Regenerate these from the EXACT run for the paper:\n")
for (nn in c(300,5000)) {
  x <- hn[hn$nn==nn & hn$method %in% c("model","robust","cluster","boot_percentile"),]
  a <- aggregate(cbind(coverage, se_ratio)~method, x, mean)
  cat(sprintf("  n=%d\n", nn)); print(a, digits=4, row.names=FALSE)
}
