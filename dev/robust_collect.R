# Collect the ARC robustness re-run and SCORE it against the exact truths.
#
# Run on the WORKSTATION, after bringing back simulation/output/robust/.
#   Rscript dev/robust_collect.R
#
# This is the scoring half of the split that dev/robust_arc.R exists to create: the cluster
# produced per-replicate estimates and no truth; every truth-dependent quantity is computed HERE.
# So re-scoring against a revised truth costs seconds and never another cluster job.
#
# It reports each cell TWICE -- against the frozen single-draw truth the original study scripts
# used, and against the exact truth that Web Tables S4 and S7 report -- because the deliverable
# is not the new number, it is the DIFFERENCE. Anything that does not move is as informative as
# anything that does.
suppressMessages({library(timereg); library(survival)})
getc <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) x else d }
ROOT <- getc("REPO_ROOT", ".")
OUT  <- getc("ROBUST_OUT", file.path(ROOT, "simulation", "output", "robust"))
A0S  <- c(0.15, 0.25, 0.40)
Z    <- qnorm(0.975)

# ONE canonical alpha key, used for every truth map and every lookup.
# format() pads to a common width on a VECTOR ("0.40") but not on a SCALAR ("0.4"), so keying
# maps with format(A0S) and reading them with format(scalar) builds a map that cannot be read
# back -- and `[[` on a missing name is a hard ERROR, not NULL, so it aborts the whole run at
# the first alpha=0.40 cell, before anything is ever saved. Three spellings coexisted across the
# two truth files and the originals ("0.40", "0.4", as.character -> "0.4").
# as.numeric() is required: this is applied to names(), which are character.
akey <- function(a) sprintf("%.2f", as.numeric(a))

# Mirrors study_grid()/STUDIES in dev/robust_arc.R. The two files drift independently, so the
# completeness check asserts against these rather than inferring them from the data it is checking.
EXPECTED_CELLS <- c(`2confounder` = 6L, misspec = 9L, tvfx = 9L)
EXPECTED_REPS  <- c(`2confounder` = 300L, misspec = 500L, tvfx = 500L)

fs <- list.files(OUT, pattern = "^rob_.*[.]rds$", full.names = TRUE)
if (!length(fs)) stop("no chunks in ", OUT)
z <- lapply(fs, readRDS)
cat(sprintf("collected %d chunk files\n", length(fs)))

# --- reassemble per-replicate rows ------------------------------------------------------------
# Column names differ by study (2confounder uses model/robust/cluster; the other two use
# se_model/se_row/se_cl). Normalise to se_model/se_row/se_cl so scoring is uniform.
norm <- c(model = "se_model", robust = "se_row", cluster = "se_cl",
          se_model = "se_model", se_row = "se_row", se_cl = "se_cl")
rows <- do.call(rbind, lapply(z, function(x) {
  e <- x$est; colnames(e) <- c("gamma", norm[x$cols[-1]])
  data.frame(study = x$study, gi = x$gi, n = x$n, alpha = x$alpha, rep = x$rep_idx,
             seed = x$seed, reps = x$reps, chunk = x$chunk, e, stringsAsFactors = FALSE)
}))

# --- completeness, BEFORE any averaging -------------------------------------------------------
cat("\n=== completeness ===\n")
bad <- FALSE
for (s in unique(rows$study)) {
  r <- rows[rows$study == s, ]
  ncell <- length(unique(r$gi))
  if (ncell != EXPECTED_CELLS[[s]]) {
    bad <- TRUE
    cat(sprintf("  !! %s: %d cells present, expected %d\n", s, ncell, EXPECTED_CELLS[[s]]))
  }
  for (g in sort(unique(r$gi))) {
    rr <- r[r$gi == g, ]
    # `want` MUST come from the reps tag each chunk carries, NOT from max(rr$rep). Inferring the
    # expected count from the data being checked means a lost trailing chunk lowers max() in
    # lockstep and the cell reports clean -- which is the dominant failure mode of a 357-task
    # array (at n=5000 a cell is 50 chunks of 10, so one walltime kill needs no operator error).
    w <- unique(rr$reps)
    if (length(w) != 1L) {
      bad <- TRUE
      cat(sprintf("  !! %s gi=%d: mixed reps tags {%s} -- stale smoke chunks?\n",
                  s, g, paste(w, collapse = ",")))
    }
    want <- w[1]
    if (want != EXPECTED_REPS[[s]]) {
      bad <- TRUE
      cat(sprintf("  !! %s gi=%d: reps tag %d, registry says %d\n",
                  s, g, want, EXPECTED_REPS[[s]]))
    }
    miss <- setdiff(seq_len(want), rr$rep)
    dup  <- sum(duplicated(rr$rep))
    if (length(miss) || dup) {
      bad <- TRUE
      cat(sprintf("  !! %s gi=%d: %d of %d replicates missing, %d duplicated\n",
                  s, g, length(miss), want, dup))
    }
  }
  cat(sprintf("  %-12s %d cells, %d replicates, %d finite gamma (%d failed)\n",
              s, ncell, nrow(r), sum(is.finite(r$gamma)), sum(!is.finite(r$gamma))))
}
if (!bad) cat("  no gaps, no duplicate replicate indices, all cells present\n")
stopifnot(!anyDuplicated(rows[, c("study", "gi", "rep")]))
# Refuse to score an incomplete collection by default. Not unconditionally fatal -- inspecting a
# run in progress is legitimate -- but it must be a deliberate act, and the result must carry the
# fact so a truncated scoring can never be mistaken for a final one.
if (bad && !nzchar(Sys.getenv("ROBUST_ALLOW_PARTIAL")))
  stop("incomplete collection; set ROBUST_ALLOW_PARTIAL=1 to score anyway")

# --- truths -----------------------------------------------------------------------------------
# EXACT for all three studies. tvfx was the last to land (dev/robustness_truth_tvfx.R): a closed
# form is indeed meaningless there -- 53% of the treated arm's k=5 hazard is truncated at
# alpha0=0.15 -- but the truncation is not an error term to be bounded, it IS the process, so
# S_a(t) = E[exp(-sum_k w_k max(h_k,0))] is evaluated directly by quadrature and the guard is
# handled exactly. There is therefore no gamma_recommended/guard-correction split for tvfx.
# Resolved through ROOT, not by bare relative path: run from another cwd (or from the shipped
# ARC bundle, which carries no truth files) and a bare path silently misses, leaving `exact`
# empty and the run completing normally with "--" everywhere -- which reads as "nothing moved",
# the exact opposite of the conclusion. Absence is now an ERROR for ALL THREE studies: with
# every truth derived, a missing file is a broken checkout, never a legitimate fall-through.
F_2CONF <- file.path(ROOT, "dev", "robustness_truth_exact.rds")
F_MISSP <- file.path(ROOT, "dev", "truth_misspec_exact.rds")
F_TVFX  <- file.path(ROOT, "dev", "truth_tvfx_exact.rds")
stopifnot(file.exists(F_2CONF), file.exists(F_MISSP), file.exists(F_TVFX))

exact <- list()
{
  e <- readRDS(F_2CONF)
  g <- if (!is.null(e$gamma_recommended)) e$gamma_recommended else e$gamma_exact
  if (is.null(e$gamma_recommended))
    cat("\n  NOTE: robustness_truth_exact.rds has no gamma_recommended; using the RAW closed form,\n",
        "  which reverts the alpha0=0.15 haz>0 guard correction (+4.1e-5). Re-run the deriver.\n")
  # The documented guard correction is +4.1e-5 at alpha0=0.15 and <=5e-10 elsewhere, so 1e-3 is
  # ~24x headroom. This exists because a partially-populated gamma_recommended once shipped with
  # element 3 left at exactly 0 (a format() scalar-vs-vector key bug in the deriver): the
  # is.null() check above cannot see that, and scoring against a truth of 0 would have produced
  # a dramatic fabricated finding.
  stopifnot(length(g) == length(e$gamma_exact), all(is.finite(g)),
            all(abs(g - e$gamma_exact) < 1e-3))
  exact[["2confounder"]] <- setNames(g, akey(e$alpha0))
}
{
  e <- readRDS(F_MISSP)
  exact[["misspec"]] <- setNames(as.numeric(unlist(e)), akey(names(e)))
}
{
  # truth_tvfx_exact.rds is a named list keyed by sprintf("%.2f", alpha0): "0.15","0.25","0.40".
  # unlist() would preserve those names, but akey(names(e)) is applied anyway so that this
  # lookup cannot depend on the deriver's spelling -- the same discipline as the misspec block.
  e <- readRDS(F_TVFX)
  exact[["tvfx"]] <- setNames(as.numeric(unlist(e)), akey(names(e)))
}

# The frozen single-draw truths the original study scripts used. Recomputed from each study's own
# true_gamma()/true_gamma2() at its hard-coded seed 2024, NOT retyped -- retyping is how a
# comparison silently becomes a comparison against a typo.
#
# CACHED, because each call is a 3e4-6e4-subject pseudo-RCT plus an aalen fit, and nine of them
# made a re-score take longer than the thing it exists to make cheap. They are deterministic
# (hard-coded seed), so caching changes nothing but the wait. Delete the cache to force a
# recompute if a study's DGM or its true_gamma() ever changes.
FROZEN_CACHE <- file.path(ROOT, "dev", "robustness_frozen_truths.rds")
frozen_cache <- if (file.exists(FROZEN_CACHE)) readRDS(FROZEN_CACHE) else list()

frozen_of <- function(study) {
  if (!is.null(frozen_cache[[study]])) return(frozen_cache[[study]])
  cat(sprintf("  computing frozen truth for %s (slow, cached afterwards)...\n", study))
  v <- frozen_compute(study)
  frozen_cache[[study]] <<- v
  saveRDS(frozen_cache, FROZEN_CACHE)
  v
}

frozen_compute <- function(study) {
  sf <- file.path(ROOT, list(`2confounder` = "dev/robustness_2confounder.R",
                             misspec = "dev/robustness_misspec.R",
                             tvfx = "dev/robustness_tvfx.R")[[study]])
  src <- readLines(sf, warn = FALSE)
  tmp <- tempfile(fileext = ".R")
  writeLines(src[1:(grep("^grid <- expand.grid", src)[1] - 1L)], tmp)
  # Private environment, NOT the global one -- the study scripts run
  # `for (f in list.files("steCI/R", ...)) source(f)`, and sourcing that globally lets the loop
  # variable clobber caller-side names. In dev/robust_arc.R that overwrote a package source file
  # before the smoke run caught it. Same hazard here; same fix.
  ENV <- new.env(parent = globalenv())
  owd <- setwd(ROOT); on.exit(setwd(owd), add = TRUE)
  suppressMessages(sys.source(tmp, envir = ENV)); setwd(owd)
  fn <- get(if (study == "2confounder") "true_gamma2" else "true_gamma", envir = ENV)
  setNames(vapply(A0S, function(a) fn(a), 0), akey(A0S))
}

# --- score ------------------------------------------------------------------------------------
cat("\n=== coverage: frozen single-draw truth vs EXACT truth ===\n")
cat("  SE ratios and Delta/V_cl are TRUTH-FREE and cannot move; they are printed once.\n\n")
res <- list()
for (s in unique(rows$study)) {
  cat(sprintf("--- %s ---\n", s))
  fr <- frozen_of(s)
  ex <- exact[[s]]
  # Fail at cell 1 with a clear message rather than at cell 5 with "subscript out of bounds".
  stopifnot(all(akey(A0S) %in% names(fr)))
  if (!is.null(ex)) stopifnot(all(akey(A0S) %in% names(ex)))
  if (is.null(ex)) cat("  (no exact truth derived for this study -- frozen column only)\n")
  cat(sprintf("%-6s %-6s %6s %9s %9s %9s   %9s %9s %9s   %8s %8s %8s %7s\n",
              "n", "alpha", "nvalid", "covM_fr", "covR_fr", "covC_fr",
              "covM_ex", "covR_ex", "covC_ex", "serM", "ser_row", "ser_cl", "D/Vcl"))
  for (g in sort(unique(rows$gi[rows$study == s]))) {
    r <- rows[rows$study == s & rows$gi == g, ]
    # EACH study's OWN filter. robustness_2confounder.R:99 keeps on is.finite(gamma) alone;
    # misspec/tvfx add is.finite(se_cl) & se_cl > 0. Imposing one uniform filter changes the
    # denominator and hence sd(gamma) and both SE ratios -- and those truth-free quantities are
    # the invariance check for the ENTIRE re-run, so breaking them for an unrelated filtering
    # reason would read as the re-run itself being unsound.
    keep <- if (s == "2confounder") is.finite(r$gamma)
            else is.finite(r$gamma) & is.finite(r$se_cl) & r$se_cl > 0
    ndrop_alt <- sum(is.finite(r$gamma)) - sum(keep)   # 0 => the two filters agree in practice
    r <- r[keep, ]
    a <- akey(r$alpha[1])
    cv <- function(tg, se) mean(abs(r$gamma - tg) <= Z * r[[se]])
    sdg <- sd(r$gamma)
    fr_c <- vapply(c("se_model", "se_row", "se_cl"), function(x) cv(fr[[a]], x), 0)
    ex_c <- if (is.null(ex)) rep(NA_real_, 3) else
            vapply(c("se_model", "se_row", "se_cl"), function(x) cv(ex[[a]], x), 0)
    # ser_model included: robustness_2confounder.R:101 publishes THREE SE ratios, so omitting
    # the model one left a third of the truth-free invariance check with no counterpart.
    ser  <- c(mean(r$se_model) / sdg, mean(r$se_row) / sdg, mean(r$se_cl) / sdg)
    dvc  <- mean(1 - (r$se_row / r$se_cl)^2)
    # ser is now (model, row, cl) -- print [2] and [3] for the row/cluster columns.
    cat(sprintf("%-6d %-6s %6d %9.3f %9.3f %9.3f   %9s %9s %9s   %8.2f %8.2f %8.2f %7.3f%s\n",
                r$n[1], a, nrow(r), fr_c[1], fr_c[2], fr_c[3],
                ifelse(is.na(ex_c[1]), "--", sprintf("%.3f", ex_c[1])),
                ifelse(is.na(ex_c[2]), "--", sprintf("%.3f", ex_c[2])),
                ifelse(is.na(ex_c[3]), "--", sprintf("%.3f", ex_c[3])),
                ser[1], ser[2], ser[3], dvc,
                if (ndrop_alt != 0) sprintf("  [filters differ by %d]", ndrop_alt) else ""))
    res[[length(res) + 1L]] <- data.frame(
      study = s, gi = g, n = r$n[1], alpha = r$alpha[1], nvalid = nrow(r),
      cov_model_frozen = fr_c[1], cov_row_frozen = fr_c[2], cov_cl_frozen = fr_c[3],
      cov_model_exact = ex_c[1], cov_row_exact = ex_c[2], cov_cl_exact = ex_c[3],
      ser_model = ser[1], ser_row = ser[2], ser_cl = ser[3],
      delta_frac = dvc, ndrop_filter_diff = ndrop_alt,
      truth_frozen = fr[[a]], truth_exact = if (is.null(ex)) NA_real_ else ex[[a]],
      mcse_cov = sqrt(0.25 / nrow(r)), row.names = NULL)
  }
  cat("\n")
}
out <- do.call(rbind, res)

cat("=== how much did the truth change actually move coverage? ===\n")
d <- out[!is.na(out$cov_cl_exact), ]
if (nrow(d)) {
  dd <- abs(c(d$cov_model_exact - d$cov_model_frozen, d$cov_row_exact - d$cov_row_frozen,
              d$cov_cl_exact - d$cov_cl_frozen))
  cat(sprintf("  max |coverage change| = %.4f over %d cell-methods\n", max(dd), length(dd)))
  cat(sprintf("  median                = %.4f\n", median(dd)))
  cat(sprintf("  typical per-cell MCSE = %.4f  <- compare against the line above\n",
              median(d$mcse_cov)))
  cat(sprintf("  cells where the change exceeds its own MCSE: %d of %d\n",
              sum(dd > rep(d$mcse_cov, 3)), length(dd)))
}
saveRDS(list(per_replicate = rows, cells = out, partial = bad),
        file.path(OUT, "robust_scored.rds"))
cat(sprintf("\nwrote %s\n", file.path(OUT, "robust_scored.rds")))
cat("per_replicate is the artifact that makes any FUTURE truth change a re-score, not a re-run.\n")
