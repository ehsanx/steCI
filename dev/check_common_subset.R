# G10/G42: interval construction fails unevenly across methods, so each method's coverage is
# computed on a slightly different replicate subset. The disclosure exists; what is missing is
# evidence the exclusion is IGNORABLE. Re-score every E2 method on the COMMON subset -- the
# replicates where all nine constructed -- and report the change.
P   <- "data/exact/results"
rep <- readRDS(file.path(P, "replicates_all.rds"))
em  <- read.csv(file.path(P, "metrics_empirical.csv"))
rep$n <- as.integer(sub("^n(\\d+)_.*", "\\1", rep$scenario))
Z <- qnorm(0.975)

tr <- em[em$estimand == "mrd", c("scenario", "tau", "truth")]
key <- function(s, t) paste(s, t)
tmap <- setNames(tr$truth, key(tr$scenario, tr$tau))

e2 <- rep[rep$estimand == "mrd" & rep$method != "point", ]
e2 <- e2[is.finite(e2$tau), ]
meths <- sort(unique(e2$method))
cat(sprintf("E2 methods: %d ; replicate-rows: %d\n", length(meths), nrow(e2)))

e2$ok  <- is.finite(e2$estimate) & is.finite(e2$se) & e2$se > 0
e2$cov <- as.integer(abs(e2$estimate - tmap[key(e2$scenario, e2$tau)]) <= Z * e2$se)

# common subset: (scenario, tau, rep) where EVERY method constructed.
# A failed method is ABSENT from the data, not present with NA, so requiring all(ok) over the
# rows that exist is trivially true and detects nothing. Require the full method COUNT.
e2$k <- paste(e2$scenario, e2$tau, e2$rep)
n_ok <- tapply(e2$ok, e2$k, sum)
common <- names(n_ok)[n_ok == length(meths)]
cat(sprintf("replicate-horizons where all %d methods constructed: %d of %d (%.2f%%)\n",
            length(meths), length(common), length(n_ok),
            100 * length(common) / length(n_ok)))
cat(sprintf("  distribution of methods-constructed per replicate-horizon: %s\n",
            paste(sprintf("%d:%d", as.integer(names(table(n_ok))), as.integer(table(n_ok))),
                  collapse = "  ")))

cat(sprintf("\n%-24s %6s %10s %10s %10s\n", "method", "n", "as-is", "common", "change"))
worst <- 0
for (m in meths) for (nn in c(300, 1000, 5000)) {
  x <- e2[e2$method == m & e2$n == nn & e2$ok, ]
  y <- x[x$k %in% common, ]
  if (!nrow(x) || !nrow(y)) next
  a <- mean(x$cov); b <- mean(y$cov)
  worst <- max(worst, abs(a - b))
  if (abs(a - b) > 0.0005)
    cat(sprintf("%-24s %6d %10.4f %10.4f %+10.4f\n", m, nn, a, b, b - a))
}
cat(sprintf("\n  worst |change| over all method x n: %.4f\n", worst))
cat("  (only cells changing by more than 0.0005 are listed above)\n")
