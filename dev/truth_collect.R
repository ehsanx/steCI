# Collect the ARC truth array and report the two things the truth check needs.
#   Rscript dev/truth_collect.R
#
# (1) VARIANCE: the averaged gamma_true per alpha.0, with its Monte-Carlo SE, against the
#     shipped single-draw value and against an E1 bias of ~0.0057.
# (2) n_hd BIAS: whether E[gamma_hat(n_hd)] is stable in n_hd. Seed-averaging cannot detect
#     this -- it converges to the expectation AT a given n_hd. If the estimate drifts across
#     the sweep, the target is n_hd-dependent and 3e4 is not enough, whatever K is.
getc <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) x else d }
ROOT  <- getc("REPO_ROOT", ".")
OUT   <- getc("TRUTH_OUT", file.path(ROOT, "simulation", "output", "truth"))
# The canonical cache lives at different paths on ARC (where only the shipped copy exists) and on the
# workstation. Try both, then give up gracefully rather than erroring --
# the "vs cached" column is informative, not load-bearing.
CANON <- Sys.getenv("TRUTH_CACHE")
if (!nzchar(CANON) || !file.exists(CANON)) {
  # dev/truth_cache_exact.rds FIRST: it holds the closed-form truths derived in
  # dev/truth_closed_form.R (checked by dev/verify_closed_form.R), and is the only
  # candidate here that is not the superseded single frozen MC draw. That draw sits
  # +0.0012 / +0.0023 / +0.0031 from exact at alpha.0 = 0.15 / 0.25 / 0.40, so a
  # 'vs cached' column computed against it measures the OLD error, not agreement
  # with the target.
  # dev/truth_cache_canonical.rds exists only inside the ARC bundle, where
  # dev/make_arc_truth_zip.py copies the canonical cache under that name; it is not in
  # this repo. The exact cache is likewise not in the bundle, so on ARC the first
  # candidate simply misses and resolution is unchanged.
  cand <- c(file.path(ROOT, "dev", "truth_cache_exact.rds"),
            file.path(ROOT, "dev", "truth_cache_canonical.rds"),
            file.path(ROOT, "data/truth_cache.rds"),
            file.path(ROOT, "simulation/output/truth_cache.rds"))
  hit  <- cand[file.exists(cand)]
  CANON <- if (length(hit)) hit[1] else ""
}
if (nzchar(CANON)) cat("[collect] canonical cache:", CANON, "\n") else
  cat("[collect] NOTE: no canonical cache found; the 'vs cached' column will read NA.\n")

fs <- list.files(OUT, pattern = "^truth_.*[.]rds$", full.names = TRUE)
if (!length(fs)) stop("no truth chunks in ", OUT)
z <- lapply(fs, readRDS)
cat(sprintf("collected %d chunk files, %d draws total\n\n", length(fs), sum(sapply(z, function(x) length(x$gamma)))))

d <- do.call(rbind, lapply(z, function(x)
  data.frame(alpha0 = x$alpha0, n_hd = x$n_hd, gamma = x$gamma)))
d <- d[is.finite(d$gamma), ]

cached <- if (nzchar(CANON) && file.exists(CANON)) {
  tc <- readRDS(CANON); sapply(names(tc), function(k) tc[[k]]$gamma_true)
} else setNames(rep(NA_real_, 3), c("0.15","0.25","0.4"))

cat("=========== (1) averaged target, by alpha.0 and n_hd ===========\n")
cat(sprintf("%-7s %8s %6s %12s %12s %12s %12s\n",
            "alpha.0", "n_hd", "K", "mean", "sd(1 draw)", "MCSE", "vs cached"))
for (a in sort(unique(d$alpha0))) for (h in sort(unique(d$n_hd))) {
  g <- d$gamma[d$alpha0 == a & d$n_hd == h]
  if (!length(g)) next
  cc <- cached[[as.character(a)]]
  cat(sprintf("%-7s %8d %6d %12.6f %12.6f %12.6f %+12.6f\n",
              format(a), h, length(g), mean(g), sd(g), sd(g)/sqrt(length(g)),
              mean(g) - cc))
}
cat("\n  E1 |bias| is ~0.005-0.0065 -- read MCSE against that.\n")

cat("\n=========== (2) is the target STABLE IN n_hd? ===========\n")
cat("If mean(gamma) drifts across n_hd beyond its MCSE, the target is n_hd-dependent and\n")
cat("seed-averaging at 3e4 converges to the WRONG number no matter how large K is.\n\n")
hs <- sort(unique(d$n_hd))
if (length(hs) < 2) {
  cat(sprintf("  only one n_hd present (%d) -- expected for the 'quick' profile.\n", hs[1]))
  cat("  The 'full' sweep is not needed for the paper's claims: the closed form\n")
  cat("  (dev/truth_closed_form.R) supplies the population gamma directly, so the fit at a\n")
  cat("  single n_hd is compared to IT rather than to other n_hd. At n_hd=3e4 the K=2000\n")
  cat("  mean agrees to |z| < 1 at every alpha.0 (dev/arc_collect_exact.R).\n")
  cat("  Read that as an UPPER BOUND on the n_hd=3e4 bias, not as proof of none: at a single\n")
  cat("  n_hd an O(1/n) bias cannot be separated from an error in the closed form.\n")
} else {
  for (a in sort(unique(d$alpha0))) {
    m <- sapply(hs, function(h) { g <- d$gamma[d$alpha0 == a & d$n_hd == h]
                                  if (length(g)) mean(g) else NA_real_ })
    s <- sapply(hs, function(h) { g <- d$gamma[d$alpha0 == a & d$n_hd == h]
                                  if (length(g)) sd(g)/sqrt(length(g)) else NA_real_ })
    cat(sprintf("  alpha.0=%-5s %s\n", format(a),
                paste(sprintf("n_hd=%d: %.6f (%.6f)", hs, m, s), collapse = "   ")))
    if (sum(is.finite(m)) >= 2) {
      i <- which(is.finite(m)); dr <- m[max(i)] - m[min(i)]
      se <- sqrt(s[max(i)]^2 + s[min(i)]^2)
      cat(sprintf("      drift %d -> %d: %+.6f  (%.1f SE)  %s\n", hs[min(i)], hs[max(i)], dr, abs(dr)/se,
                  if (abs(dr) > 3*se) "<-- n_hd-DEPENDENT, use the largest n_hd" else "stable"))
    }
  }
}

cat("\n=========== gamma_true is EXACT -- do NOT cache the numbers below ===========\n")
cat("  The E1 target has a closed form derived from the DGM (dev/truth_closed_form.R,\n")
cat("  checked by dev/verify_closed_form.R). The exact values live in the repo at\n")
cat("  dev/truth_cache_exact.rds:\n")
cat("      alpha.0 = 0.15  ->  -0.06183507\n")
cat("      alpha.0 = 0.25  ->  -0.05965809\n")
cat("      alpha.0 = 0.40  ->  -0.05656881\n")
cat("  The means below are Monte-Carlo estimates carrying MCSE; adopting them as the\n")
cat("  target would replace an exact number with a noisier one. They are useful only\n")
cat("  as a comparison against the exact values, which dev/arc_collect_exact.R does\n")
cat("  properly. (That script and the exact cache live in the repo, not in the ARC\n")
cat("  zip -- bring the draws home and run it there.)\n")
best <- max(d$n_hd)
cat(sprintf("MC mean at the largest n_hd present (%d), for COMPARISON only:\n", best))
for (a in sort(unique(d$alpha0))) {
  g <- d$gamma[d$alpha0 == a & d$n_hd == best]
  if (length(g)) cat(sprintf("  alpha.0 = %-5s MC mean = %.6f   (MCSE %.6f, K = %d)\n",
                             format(a), mean(g), sd(g)/sqrt(length(g)), length(g)))
}
saveRDS(d, file.path(OUT, "truth_draws_all.rds"))
cat("\nwrote", file.path(OUT, "truth_draws_all.rds"), "\n")
