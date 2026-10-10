# Empirical population limit of the POOLED weighted-Aalen STE estimator, on ARC.
#
# WHY THIS EXISTS. Web-review Finding 2 (round 2) asks whether the residual long-horizon E2
# point-estimate bias is a NON-VANISHING pooled-origin offset (the pooling convention B(.) is
# common across trial origins, which the DGM violates) or a finite-sample effect that shrinks
# with n. Three sample sizes (n=300/1000/5000) cannot separate a vanishing bias from a plateau.
# This job extends the trend to LARGE n at high precision, so the empirical plim of the SHIPPED
# estimator can be read off and cross-checked against the closed-form pooled plim
# (dev/plim_closed_form.R). If the two agree, the pooled offset is pinned exactly.
#
# THE DESIGN POINT: THE TRUTH IS AN INPUT TO SCORING, NOT TO SIMULATION (repo standing rule 13).
# This script computes NO truth and reads NO truth. It emits one row per replicate --
# (gamma_hat, mrd_hat[1..5]) -- and stops. Scoring (bias vs the exact closed-form truths, SEs,
# the vanish-vs-plateau read, and the closed-form cross-check) happens at HOME in
# dev/plim_collect.R. Any future truth revision is then a seconds-long re-score, never a re-run.
#
# THE ESTIMATOR IS NOT REIMPLEMENTED. It calls the shipped ste_estimate() on dat_sim() output,
# byte-identical to the code behind the published tables, so a disagreement with the closed form
# is a derivation issue, not a transcription slip.
#
# One array task = one (cell, n, replicate-chunk) -> one RDS. Checkpointed on file existence.
# Env: PLIM_CELLS (colon-separated cell ids, e.g. "1:2:3:4"), PLIM_NS (colon-separated n),
#      PLIM_REPS_OVERRIDE (smoke only), PLIM_OUT, REPO_ROOT.
suppressMessages({library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang)})

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
getc <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) x else d }
ROOT <- getc("REPO_ROOT", ".")
OUT  <- getc("PLIM_OUT", file.path(ROOT, "simulation", "output", "plim"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

owd <- setwd(ROOT); on.exit(setwd(owd), add = TRUE)   # steCI/R and scenarios.R by relative path
for (f in list.files(file.path(ROOT, "steCI", "R"), pattern = "[.]R$", full.names = TRUE)) source(f)
# scenarios.R gives DGM_FIXED and T_HOR; the two repo copies are byte-identical, so take whichever
# ships (the ARC zip carries simulation/scenarios.R). Fail loudly if neither is present.
scen <- Filter(file.exists, file.path(ROOT, c(
  "simulation/scenarios.R",
  "simulation/scenarios.R")))
stopifnot(length(scen) >= 1L); source(scen[1])
setwd(owd)
stopifnot(is.function(dat_sim), is.function(ste_estimate), is.list(DGM_FIXED), length(T_HOR) == 5L)

# --- CELLS ------------------------------------------------------------------------------------
# The strong-confounding low-event-rate corner is where the observed E2 bias lives (worst cell
# n5000_evlow_cfstrong_pvhigh, +0.0248 at tau=5); the weak-confounding twin is the ~zero-bias
# control; the moderate cell fills the confounding axis; the high-event-rate strong cell probes
# the event-rate dependence. alpha.0/gamma.L/gamma.0 map to event/confounding/prevalence exactly
# as in scenarios.R (event low/med/high = alpha.0 0.15/0.25/0.40; conf weak/mod/strong = gamma.L
# 0.5/1.5/3.0; prev high/med/low = gamma.0 0/-1/-2). The E2 truth depends on alpha.0 ONLY.
CELLS <- data.frame(
  cell_id = 1:4,
  label   = c("evlow_cfstrong_pvhigh", "evlow_cfweak_pvhigh",
              "evlow_cfmod_pvhigh",    "evhigh_cfstrong_pvhigh"),
  alpha.0 = c(0.15, 0.15, 0.15, 0.40),
  gamma.L = c(3.0,  0.5,  1.5,  3.0),
  gamma.0 = c(0.0,  0.0,  0.0,  0.0),
  stringsAsFactors = FALSE)

NS   <- as.integer(strsplit(getc("PLIM_NS", "10000"), "[:,]")[[1]])
NS   <- NS[!is.na(NS)]
# n IS CAPPED AT 1e4 BY DESIGN (2026-07-23). Two things forced this:
#   (1) timereg::aalen is empirically ~O(n^2) here -- MEASURED on ARC (smoke job 12342697): 85-122
#       s/fit at n=1e4 and 985 s/fit at n=3e4 -- so n=1e5 is ~3 h/fit (its four smoke tasks were all
#       cancelled at the walltime) and even n=3e4 is 16 min/fit.
#   (2) MORE IMPORTANTLY, the n=3e4 arm carries NO RESOLVABLE SIGNAL. The finite-n bias moves only
#       ~0.002 between n=5000 and n=3e4 (the curve has flattened toward the closed-form limit
#       +0.007), while emp_sd at tau=5 is ~0.072 (heavy-tailed weights), so 80 reps give SE ~0.008
#       and "declining toward 0.007" vs "stays at 0.025" is z~0.2. A power analysis
#       shows NO feasible replicate count rescues n>=3e4 (3-sigma would need ~1e4 reps/cell, and the
#       replicate mean is itself heavy-tailed so SD/sqrt(R) is not even a valid SE).
# So this array does ONE useful thing: an INDEPENDENT-CODE cross-check that the SHIPPED ste_estimate
# (real timereg::aalen, no Rao-Blackwellisation) gives a clearly POSITIVE finite-n bias at n=1e4
# (~0.024 +/- 0.006), consistent with the closed form's positive asymptote -- closing the "is your
# closed form the same object as your package?" gap. The n -> infinity value stays with
# dev/plim_closed_form.R. To run n=3e4 anyway (not recommended), set PLIM_NS=10000:30000.
reps_for <- function(n) if (n >= 30000L) 80L else 200L        # ~0.006 SE at n=1e4/200
chunk_reps <- function(n) if (n >= 30000L) 1L else 5L         # 5 fits/task (~10 min) at n=1e4

# --- flat task list ---------------------------------------------------------------------------
want_cells <- as.integer(strsplit(getc("PLIM_CELLS", "1:2:3:4"), "[:,]")[[1]])
want_cells <- want_cells[!is.na(want_cells)]
stopifnot(all(want_cells %in% CELLS$cell_id))
REPS_OVR <- geti("PLIM_REPS_OVERRIDE", NA_integer_)     # smoke only; NEVER set for production

tasks <- do.call(rbind, lapply(want_cells, function(cid) {
  do.call(rbind, lapply(seq_along(NS), function(ni) {
    n <- NS[ni]; R <- if (is.na(REPS_OVR)) reps_for(n) else REPS_OVR
    cr <- min(chunk_reps(n), R)
    data.frame(cell_id = cid, n_idx = ni, n = n, reps = R,
               chunk = seq_len(ceiling(R / cr)), chunk_reps = cr, stringsAsFactors = FALSE)
  }))
}))

tid <- geti("SLURM_ARRAY_TASK_ID", 1L)
cat(sprintf("[plim] task %d of %d total tasks (cells: %s ; n: %s)\n",
            tid, nrow(tasks), paste(want_cells, collapse = ","), paste(NS, collapse = ",")))
if (tid > nrow(tasks)) { cat("[plim] beyond grid -- nothing to do\n"); quit(save = "no") }
tk <- tasks[tid, ]
cell <- CELLS[CELLS$cell_id == tk$cell_id, ]

# Cell, n, reps and chunk all in the filename: a smoke run (R-tag) cannot poison production, and
# a changed grid cannot silently overwrite a previous chunk.
outfile <- file.path(OUT, sprintf("plim_c%d_%s_n%d_R%d_chunk%03d.rds",
                            tk$cell_id, cell$label, tk$n, tk$reps, tk$chunk))
if (file.exists(outfile)) { cat("[plim] already complete:", basename(outfile), "\n"); quit(save = "no") }
stopifnot(is.character(outfile), length(outfile) == 1L,
          grepl("^plim_.*[.]rds$", basename(outfile)),
          identical(normalizePath(dirname(outfile), mustWork = TRUE),
                    normalizePath(OUT, mustWork = TRUE)))

idx <- ((tk$chunk - 1L) * tk$chunk_reps) + seq_len(tk$chunk_reps)
idx <- idx[idx <= tk$reps]
stopifnot(length(idx) >= 1L)
cat(sprintf("[plim] cell %d (%s) n=%d | chunk %d, replicates %d..%d\n",
            tk$cell_id, cell$label, tk$n, tk$chunk, min(idx), max(idx)))

# est columns: gamma (E1) + mrd at the five horizons.
est  <- matrix(NA_real_, length(idx), 6, dimnames = list(NULL, c("gamma", paste0("mrd", T_HOR))))
fail <- character(length(idx))
t0   <- Sys.time()
for (i in seq_along(idx)) {
  sd_i <- 100000L * tk$cell_id + 1000L * tk$n_idx + idx[i]   # collision-free per (cell, n, rep)
  d <- tryCatch(
    dat_sim(n = tk$n, n.visit = DGM_FIXED$n.visit, gamma.0 = cell$gamma.0, gamma.L = cell$gamma.L,
            alpha.0 = cell$alpha.0, alpha.A = DGM_FIXED$alpha.A, alpha.L = DGM_FIXED$alpha.L,
            alpha.U = DGM_FIXED$alpha.U, seed = sd_i),
    error = function(e) { fail[i] <<- paste("dat_sim:", conditionMessage(e)); NULL })
  if (is.null(d)) next
  e <- tryCatch(ste_estimate(d, n.visit = DGM_FIXED$n.visit, t.hor = T_HOR,
                             weight = "stab", estimand = "both"),
                error = function(e) { fail[i] <<- paste("ste_estimate:", conditionMessage(e)); NULL })
  if (!is.null(e)) est[i, ] <- c(e$gamma, e$mrd)
}
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

peak_rss_mb <- function() {
  if (!file.exists("/proc/self/status")) return(NA_real_)
  l <- grep("^VmHWM:", readLines("/proc/self/status", warn = FALSE), value = TRUE)
  if (!length(l)) return(NA_real_)
  as.numeric(sub("\\D*(\\d+).*", "\\1", l[1])) / 1024
}

# ATOMIC WRITE (see robust_arc.R): checkpointed on file existence, so a mid-write kill must not
# leave a truncated .rds that the next submission skips as complete. Temp in the SAME dir + rename.
tmp_out <- paste0(outfile, ".part", Sys.getpid())
saveRDS(list(cell_id = tk$cell_id, label = cell$label, alpha.0 = cell$alpha.0,
             gamma.L = cell$gamma.L, gamma.0 = cell$gamma.0, n = tk$n, reps = tk$reps,
             chunk = tk$chunk, rep_idx = idx, seed = 100000L * tk$cell_id + 1000L * tk$n_idx + idx,
             t.hor = T_HOR, cols = colnames(est), est = est, fail = fail,
             n_ok = sum(is.finite(est[, 1])), secs = secs, rss_mb = peak_rss_mb()), tmp_out)
if (!file.rename(tmp_out, outfile)) { unlink(tmp_out); stop("could not rename ", tmp_out, " -> ", outfile) }

cat(sprintf("[plim] wrote %s | %d/%d ok | %.1fs (%.1fs/fit) | peak %s MB\n", basename(outfile),
            sum(is.finite(est[, 1])), length(idx), secs, secs / max(1L, length(idx)),
            ifelse(is.na(peak_rss_mb()), "n/a", sprintf("%.0f", peak_rss_mb()))))
if (any(nzchar(fail))) {
  cat(sprintf("[plim] WARNING: %d of %d replicates failed\n", sum(nzchar(fail)), length(idx)))
  for (m in unique(fail[nzchar(fail)])) cat("   ", m, "\n")
}
