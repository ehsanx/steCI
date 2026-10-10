# Score the pooled-plim ARC run AT HOME (the cluster computes no truth; standing rule 13).
#
# Reads the per-replicate chunks in simulation/output/plim/, averages the SHIPPED estimator's
# gamma/mrd per (cell, n), applies the EXACT closed-form truths, and reports:
#   (1) the tau=5 E2 bias vs n for each cell, stitched onto the canonical n=300/1000/5000 trend,
#       so a vanishing bias can be told from a plateau;
#   (2) the empirical plim (mean at the largest n) with its Monte-Carlo SE;
#   (3) a cross-check against the closed-form pooled plim (dev/plim_closed_form.rds), if present.
#
#   Rscript dev/plim_collect.R
suppressWarnings(suppressMessages({library(stats)}))

getc <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) x else d }
ROOT   <- getc("REPO_ROOT", ".")
PLIMD  <- getc("PLIM_OUT", file.path(ROOT, "simulation", "output", "plim"))
CANON  <- file.path(ROOT, "data", "exact", "results", "metrics_empirical.csv")
TCACHE <- file.path(ROOT, "dev", "truth_cache_exact.rds")

files <- list.files(PLIMD, pattern = "^plim_.*\\.rds$", full.names = TRUE)
files <- files[!grepl("_R2_", files)]                     # drop smoke chunks
if (!length(files)) stop("no plim_*.rds in ", PLIMD, " (run the ARC job, then tar the output home)")
cat(sprintf("[plim] %d production chunks in %s\n", length(files), PLIMD))

# --- pool chunks into one per-replicate frame ------------------------------------------------
rows <- lapply(files, function(f) {
  x <- readRDS(f); e <- x$est
  ok <- is.finite(e[, "gamma"])                           # a valid replicate has a point estimate
  if (!any(ok)) return(NULL)
  data.frame(cell_id = x$cell_id, label = x$label, alpha.0 = x$alpha.0,
             gamma.L = x$gamma.L, n = x$n, e[ok, , drop = FALSE],
             row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
})
d <- do.call(rbind, Filter(Negate(is.null), rows))
thor <- as.integer(sub("mrd", "", grep("^mrd", colnames(d), value = TRUE)))
mcols <- paste0("mrd", thor)

# --- exact truths (E2 mrd depends on alpha.0 only) -------------------------------------------
tv <- readRDS(TCACHE)
truth_mrd <- function(a0, tau) tv[[as.character(a0)]]$mrd[match(tau, tv[[as.character(a0)]]$t.hor)]

# --- aggregate per (cell, n): mean, MC SE, bias ----------------------------------------------
key <- paste(d$cell_id, d$n, sep = "|")
agg <- do.call(rbind, lapply(split(d, key), function(g) {
  a0 <- g$alpha.0[1]; nrep <- nrow(g)
  out <- data.frame(cell_id = g$cell_id[1], label = g$label[1], alpha.0 = a0,
                    gamma.L = g$gamma.L[1], n = g$n[1], n_rep = nrep, stringsAsFactors = FALSE)
  for (tau in thor) {
    est <- g[[paste0("mrd", tau)]]; est <- est[is.finite(est)]
    m <- mean(est); se <- sd(est) / sqrt(length(est)); tr <- truth_mrd(a0, tau)
    out[[sprintf("bias_t%d", tau)]]   <- m - tr
    out[[sprintf("se_t%d", tau)]]     <- se
    out[[sprintf("mean_t%d", tau)]]   <- m
  }
  out
}))
agg <- agg[order(agg$cell_id, agg$n), ]

# --- canonical n=300/1000/5000 tau=5 bias for the same cells, to complete the trend ----------
canon_t5 <- NULL
if (file.exists(CANON)) {
  me <- read.csv(CANON, stringsAsFactors = FALSE)
  me5 <- me[me$estimand == "mrd" & me$tau == 5, ]
  canon_t5 <- do.call(rbind, lapply(unique(agg$label), function(lab) {
    do.call(rbind, lapply(c(300L, 1000L, 5000L), function(nn) {
      r <- me5[me5$scenario == sprintf("n%d_%s", nn, lab), ]
      if (!nrow(r)) return(NULL)
      data.frame(label = lab, n = nn, bias_t5 = r$bias[1],
                 se_t5 = if ("emp_sd" %in% names(r) && "n_sim" %in% names(r))
                           r$emp_sd[1] / sqrt(r$n_sim[1]) else NA_real_,
                 stringsAsFactors = FALSE)
    }))
  }))
}

# --- report: tau=5 bias vs n (the vanish-vs-plateau read) ------------------------------------
cat("\n================  E2 tau=5 bias vs n  (canonical + ARC large-n)  ================\n")
for (lab in unique(agg$label)) {
  cat(sprintf("\n%s  (alpha.0=%.2f, gamma.L=%.1f):\n", lab,
              agg$alpha.0[agg$label == lab][1], agg$gamma.L[agg$label == lab][1]))
  tr <- rbind(
    if (!is.null(canon_t5)) {
      c5 <- canon_t5[canon_t5$label == lab, ]
      if (nrow(c5)) data.frame(n = c5$n, bias = c5$bias_t5, se = c5$se_t5, src = "canon") else NULL
    },
    { a <- agg[agg$label == lab, ]; data.frame(n = a$n, bias = a$bias_t5, se = a$se_t5, src = "arc") })
  tr <- tr[order(tr$n), ]
  for (i in seq_len(nrow(tr)))
    cat(sprintf("   n=%-7d  bias=%+.5f  +/- %.5f  (%s)%s\n", tr$n[i], tr$bias[i], tr$se[i], tr$src[i],
                if (is.finite(tr$se[i]) && abs(tr$bias[i]) > 2 * tr$se[i]) "  [>2 SE from 0]" else ""))
}

# --- cross-check against the closed-form pooled plim -----------------------------------------
CF <- file.path(ROOT, "dev", "plim_closed_form.rds")
if (file.exists(CF)) {
  cf <- readRDS(CF)   # expected: data.frame(label, tau, mrd_pooled_plim [, mrd_origin_plim, truth])
  cat("\n================  cross-check: ARC empirical (largest n) vs closed-form plim  ============\n")
  for (lab in unique(agg$label)) {
    a <- agg[agg$label == lab, ]; a <- a[a$n == max(a$n), ]
    cc <- cf[cf$label == lab, ]
    if (!nrow(cc)) next
    cat(sprintf("\n%s  (empirical at n=%d, %d reps):\n", lab, a$n, a$n_rep))
    # Current dev/plim_closed_form.rds schema (RB version): plim_pooled + se_pooled. (Older objects
    # used mrd_pooled_plim/se_jack -- fall back to those so this works on either.)
    pcol <- if ("plim_pooled" %in% names(cc)) "plim_pooled" else "mrd_pooled_plim"
    scol <- if ("se_pooled"   %in% names(cc)) "se_pooled"   else "se_jack"
    for (tau in thor) {
      emp <- a[[sprintf("mean_t%d", tau)]]; se <- a[[sprintf("se_t%d", tau)]]
      pl  <- cc[[pcol]][cc$tau == tau]
      if (!length(pl)) next
      # The plim is NOT exact and its stored jackknife SE UNDERSTATES the truth at long horizons (it
      # is a within-stream measure; independent seed streams disagree ~3x at tau=5 -- see
      # dev/plim_seedcheck.R and the plim_closed_form.R header). So this z is only meaningful at
      # tau<=2; at tau>=3 treat a "MISMATCH" flag as inconclusive, not as a real disagreement.
      se_pl <- if (!is.null(cc[[scol]])) cc[[scol]][cc$tau == tau] else 0
      se_tot <- sqrt(se^2 + se_pl^2)
      z <- (emp - pl) / se_tot
      cat(sprintf("   tau=%d  empirical=%+.5f +/- %.5f   closed-form=%+.5f +/- %.5f%s   z=%+.2f%s\n",
                  tau, emp, se, pl, se_pl, if (tau >= 3) " [SE understated]" else "",
                  z, if (abs(z) > 3 && tau <= 2) "  <-- MISMATCH" else ""))
    }
  }
} else {
  cat("\n[plim] dev/plim_closed_form.rds not found -- run dev/plim_closed_form.R for the cross-check.\n")
}

saveRDS(list(agg = agg, canon_t5 = canon_t5, n_files = length(files)),
        file.path(ROOT, "dev", "plim_results.rds"))
cat("\n[plim] wrote dev/plim_results.rds\n")
