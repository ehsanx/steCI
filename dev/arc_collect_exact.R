# Collect the completed ARC truth array (K=2000, n_hd=3e4, 189 chunks) and test it
# against the EXACT closed-form truths.
#
# The run's purpose changed after M8 dissolved. It no longer ESTIMATES the target -- the
# closed form gives that exactly. It now measures whether the timereg::aalen const() fit
# used by true_values() carries finite-sample bias at n_hd = 3e4. That is a question the
# closed form cannot answer and seed-averaging alone cannot detect.
#
#   Rscript dev/arc_collect_exact.R
OUT <- "simulation/steCI_truth_arc/simulation/output/truth"
fs  <- list.files(OUT, pattern = "^truth_.*chunk[0-9]+[.]rds$", full.names = TRUE)
stopifnot(length(fs) == 189L)
z <- lapply(fs, readRDS)

# Structural checks before any averaging: a silently short or duplicated array would
# otherwise average to a plausible-looking wrong number.
stopifnot(all(vapply(z, function(x) x$K, numeric(1)) == 2000))
stopifnot(all(vapply(z, function(x) x$n_hd, numeric(1)) == 30000))
stopifnot(all(vapply(z, function(x) x$n_died, numeric(1)) == 0))
sd_all <- unlist(lapply(z, function(x) x$seed))
stopifnot(!anyDuplicated(sd_all))          # globally unique seeds, as designed

d <- do.call(rbind, lapply(z, function(x)
       data.frame(alpha0 = x$alpha0, chunk = x$chunk, seed = x$seed, gamma = x$gamma)))
cat(sprintf("chunks %d | draws %d | finite %d | unique seeds %d\n",
            length(fs), nrow(d), sum(is.finite(d$gamma)), length(unique(d$seed))))
for (a in sort(unique(d$alpha0)))
  stopifnot(sum(d$alpha0 == a) == 2000L)   # exactly K per arm, no chunk lost or doubled

ex  <- readRDS("dev/truth_cache_exact.rds")
can <- readRDS("simulation/output/truth_cache.rds")

cat("\n=== does the aalen fit at n_hd=3e4 hit the exact truth? ===\n")
cat(sprintf("%-6s %6s %12s %11s %11s %14s %+11s %7s\n",
            "alpha0", "K", "mean(ARC)", "sd(1 draw)", "MCSE", "exact", "ARC-exact", "z"))
res <- list()
for (a in sort(unique(d$alpha0))) {
  g <- d$gamma[d$alpha0 == a]; m <- mean(g); s <- sd(g); se <- s / sqrt(length(g))
  e <- ex[[as.character(a)]]$gamma_true
  res[[as.character(a)]] <- list(mean = m, sd = s, mcse = se, exact = e, z = (m - e) / se)
  cat(sprintf("%-6s %6d %12.6f %11.6f %11.6f %14.8f %+11.6f %7.2f\n",
              format(a), length(g), m, s, se, e, m - e, (m - e) / se))
}
cat("\n=== the BOUND -- report this, never 'unbiased' ===\n")
cat("  |z| < 1 accepts a null; it does not prove one. The defensible quantity is the larger\n")
cat("  95% CI endpoint |delta| + 1.96*MCSE -- NOT 2*MCSE, which ignores that delta != 0 and\n")
cat("  understates the binding arm (alpha.0 = 0.4) by 46%.\n")
for (a in sort(unique(d$alpha0))) {
  r <- res[[as.character(a)]]; b <- r$mean - r$exact; h <- 1.96 * r$mcse
  cat(sprintf("  alpha.0=%-5s CI [%+.2e, %+.2e]  ->  |bias| < %.2e\n",
              format(a), b - h, b + h, abs(b) + h))
}
worst <- max(vapply(res, function(r) abs(r$mean - r$exact) + 1.96 * r$mcse, 0))
cat(sprintf("\n  WORST ARM: |bias| < %.5f  (%.1f%% of the +0.007 E1 bias it measures)\n",
            worst, 100 * worst / 0.007))
X <- sum(vapply(res, function(r) r$z^2, 0))
cat(sprintf("  OMNIBUS:   sum z^2 = %.3f on 3 df, p = %.3f (seed blocks disjoint => independent)\n",
            X, pchisq(X, 3, lower.tail = FALSE)))
cat("\n  NOTE this is a COMPOUND null: it tests (closed form correct) AND (aalen unbiased at\n")
cat("  n_hd=3e4) jointly, and at a single n_hd the two are not separable. Lead with the\n")
cat("  CONFIRMATION reading -- the paper uses the closed form, not the aalen fit.\n")

cat("\n=== what the single frozen draw cost (M8's original complaint) ===\n")
for (a in sort(unique(d$alpha0))) {
  cv <- can[[as.character(a)]]$gamma_true; e <- ex[[as.character(a)]]$gamma_true
  cat(sprintf("  alpha.0=%-5s shipped %.8f  exact %.8f  displacement %+.6f  (%.1f single-draw SD)\n",
              format(a), cv, e, cv - e, abs(cv - e) / res[[as.character(a)]]$sd))
}

saveRDS(list(draws = d, summary = res), file.path("dev", "arc_truth_results.rds"))
cat("\nwrote dev/arc_truth_results.rds\n")
