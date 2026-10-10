# Score the max.clust sensitivity AT HOME (the cluster computes no truth; standing rule 13).
#
# Reads the per-replicate (gamma, se_default, se_pinned, nclust) chunks, applies the EXACT E1
# truths (dev/truth_cache_exact.rds), and reports:
#   (0) a BINNING SANITY check -- did timereg's default actually bin (max.clust < nclust)? and how
#       many replicates show an exactly-zero contrast? A large exact-zero fraction means the
#       default-vs-NULL comparison did not fire and the numbers below are vacuous.
#   (1) the relative M3 SE shift se_pinned/se_default - 1, LED BY MAGNITUDE (mean, mean|shift|,
#       per-cell), with a t-test reported only as a detectability check (at ~2700 reps it will
#       flag even a practically-negligible shift; on the SE scale Jensen makes the shift positive
#       in expectation, so a small significant positive value is EXPECTED and immaterial);
#   (2) the induced E1 coverage change, judged against the PAIRED (McNemar) SE -- both arms use the
#       identical reps, the pin-invariant gamma and the identical truth, so only discordant reps
#       move it; the single-arm binomial SE is the wrong, much looser yardstick.
#
#   Rscript dev/maxclust_collect.R
suppressWarnings(suppressMessages({library(stats); library(timereg)}))
getc <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) x else d }
ROOT <- getc("REPO_ROOT", ".")
DIR  <- getc("MC_OUT", file.path(ROOT, "simulation", "output", "maxclust"))
TC   <- file.path(ROOT, "dev", "truth_cache_exact.rds")

files <- list.files(DIR, pattern = "^mc_.*\\.rds$", full.names = TRUE)
files <- files[!grepl("_R2_", basename(files))]        # drop smoke chunks (match basename only)
if (!length(files)) stop("no production mc_*.rds in ", DIR, " (run the ARC job, then tar it home)")
cat(sprintf("[mc] %d production chunks in %s\n", length(files), DIR))

n_fail <- 0L
rows <- lapply(files, function(f) {
  x <- readRDS(f)
  n_fail <<- n_fail + sum(nzchar(x$fail))
  ok <- is.finite(x$se_default) & is.finite(x$se_pinned) & !nzchar(x$fail)
  if (!any(ok)) return(NULL)
  data.frame(scenario = x$scenario, alpha.0 = x$alpha.0, nclust = x$nclust[ok], gamma = x$gamma[ok],
             se_default = x$se_default[ok], se_pinned = x$se_pinned[ok],
             row.names = NULL, stringsAsFactors = FALSE)
})
d <- do.call(rbind, Filter(Negate(is.null), rows))
if (is.null(d) || !nrow(d)) stop("no valid replicate comparisons found in ", DIR)
d$rel <- d$se_pinned / d$se_default - 1

# --- (completeness) -------------------------------------------------------------------------
cells <- sort(unique(d$scenario))
cat(sprintf("\n%d valid replicate comparisons across %d n=5000 cells (%d replicate-level failures skipped)\n",
            nrow(d), length(cells), n_fail))
if (length(cells) != 27L)
  cat(sprintf("  !! expected 27 n=5000 cells, found %d -- run is INCOMPLETE; numbers are a subset.\n",
              length(cells)))

# --- (0) binning sanity: did the default actually bin, and is the contrast non-degenerate? ---
mc_def <- formals(timereg::aalen)$max.clust
frac0  <- mean(d$rel == 0)
cat("\n================  binning sanity  ================\n")
cat(sprintf("  timereg %s ; resolved default max.clust = %s ; cluster count (nclust) range %d-%d\n",
            as.character(packageVersion("timereg")), deparse(mc_def), min(d$nclust), max(d$nclust)))
cat(sprintf("  binning expected (default < nclust)? %s ; replicates with EXACTLY zero shift: %d/%d (%.0f%%)\n",
            is.numeric(mc_def) && mc_def < min(d$nclust), sum(d$rel == 0), nrow(d), 100 * frac0))
if (frac0 > 0.5)
  cat("  !! WARNING: over half the replicates show a bit-identical se_pinned == se_default. The\n",
      "     default-vs-NULL contrast is largely ABSENT (binning did not fire, or the pin was not\n",
      "     applied). Do NOT read the shift/coverage below as a confirmation -- investigate first.\n", sep = "")

# --- (1) relative M3 SE shift: magnitude first ----------------------------------------------
m <- mean(d$rel); s <- sd(d$rel); se <- s / sqrt(nrow(d)); tt <- m / se
cat("\n================  relative M3 SE shift: se_pinned / se_default - 1  ================\n")
cat(sprintf("  MAGNITUDE:  mean %+.4f%%   mean|shift| %.4f%%   sd %.3f%%   fraction positive %.2f\n",
            100 * m, 100 * mean(abs(d$rel)), 100 * s, mean(d$rel > 0)))
cat(sprintf("  detectability (NOT relevance): t = %.2f (df=%d), p = %.3g, SE of mean %.4f%%.\n",
            tt, nrow(d) - 1, 2 * pt(-abs(tt), nrow(d) - 1), 100 * se))
cat("  At this N the t-test flags even a negligible shift; on the SE scale the shift is positive in\n")
cat("  expectation (Jensen), so a small significant positive value is expected and immaterial.\n")

cat("\n== per-cell mean relative shift (%) ==\n")
ag <- do.call(rbind, lapply(split(d, d$scenario), function(g)
  data.frame(scenario = g$scenario[1], mean_pct = 100 * mean(g$rel),
             sd_pct = 100 * sd(g$rel), n = nrow(g), stringsAsFactors = FALSE)))
ag <- ag[order(ag$mean_pct), ]
print(ag, row.names = FALSE, digits = 3)
# Cell-treated-as-unit summary (equal-weight over cells), a cleaner overall than the pooled SE.
cm <- ag$mean_pct
cat(sprintf("\n  cell-level (equal-weight): mean of the %d cell means %+.4f%%, between-cell SE %.4f%%\n",
            length(cm), mean(cm), sd(cm) / sqrt(length(cm))))

# --- (2) induced E1 coverage change: PAIRED (McNemar) SE ------------------------------------
if (file.exists(TC)) {
  tv <- readRDS(TC)
  d$truth <- vapply(as.character(d$alpha.0), function(a) tv[[a]]$gamma_true, numeric(1))
  z <- qnorm(0.975)
  cov_d <- abs(d$gamma - d$truth) <= z * d$se_default   # covered under the published (default) SE
  cov_p <- abs(d$gamma - d$truth) <= z * d$se_pinned    # covered under max.clust = NULL
  di <- as.integer(cov_p) - as.integer(cov_d)           # +1 gained, -1 lost, 0 concordant
  diff <- mean(di); se_diff <- sd(di) / sqrt(length(di))
  gained <- sum(di > 0); lost <- sum(di < 0)
  cat("\n================  induced E1 coverage change at n=5000  ================\n")
  cat(sprintf("  coverage(default / published path) = %.4f\n", mean(cov_d)))
  cat(sprintf("  coverage(max.clust = NULL)         = %.4f\n", mean(cov_p)))
  cat(sprintf("  paired difference (NULL - default) = %+.4f   PAIRED SE %.4f   (discordant: +%d / -%d of %d)\n",
              diff, se_diff, gained, lost, nrow(d)))
  cat("  Judge the difference against the PAIRED SE (only discordant reps move it), not a single-arm\n")
  cat("  binomial SE. A difference within a few paired SEs leaves the published entry effectively unmoved.\n")
} else {
  cat("\n[mc] dev/truth_cache_exact.rds not found -- coverage change skipped (run at home).\n")
  diff <- se_diff <- NA_real_
}

saveRDS(list(rel = d$rel, per_cell = ag, n_rep = nrow(d), n_cells = length(cells),
             frac_exact_zero = frac0, mean_shift = m, mean_abs_shift = mean(abs(d$rel)),
             cov_diff = diff, cov_diff_se = se_diff),
        file.path(ROOT, "dev", "maxclust_results.rds"))
cat("\n[mc] wrote dev/maxclust_results.rds\n")
