# max.clust sensitivity on ARC. Does pinning timereg's max.clust = NULL change the M3
# cluster-robust SE of the constant hazard difference (E1)?
#
# WHY THIS EXISTS. The published M3 cluster-robust SE is timereg::aalen(..., clusters = id) at
# timereg's DEFAULT max.clust = 1000. When the number of subject clusters exceeds 1000, timereg
# quantile-bins them into 1000 groups, so at n = 5000 (5000 clusters > 1000) distinct subjects
# share bins; at n = 300 and n = 1000 the clusters are exact (1000 is not < 300/1000), so pinning
# would change nothing there. The binned "meat" is UNBIASED for the exact subject-level meat
# (subjects are independent, so the omitted cross-cluster term is mean-zero), so the cost is extra
# noise, not bias -- but the claim rests on a 20-replicate check (mean +0.37%, sd 1.63%, p = 0.33).
# This job tightens that bound: it re-fits the cluster SE with and without the pin over many
# replicates at every n = 5000 cell, so the mean shift and its induced coverage change are pinned
# down. n < 5000 is exact by construction and is not run.
#
# THE ESTIMATOR IS NOT REIMPLEMENTED and NO TRUTH IS COMPUTED OR READ. The
# default fit is byte-identical to ci_analytic.R's cluster fit (verified locally to reproduce the
# published SE); max.clust changes ONLY the robust variance, never gamma (asserted below). Each
# replicate emits (gamma, se_default, se_pinned); dev/maxclust_collect.R applies the exact E1
# truths AT HOME and reports the relative SE shift, its t-test, and the induced coverage change.
#
# One array task = one (cell, replicate-chunk) -> one RDS. Checkpointed on file existence.
# Env: MC_CELLS (colon-separated 1-based indices into the 27 n=5000 cells; empty = all),
#      MC_REPS (per cell, default 100), MC_CHUNK (reps per task, default 10),
#      MC_REPS_OVERRIDE (smoke only; NEVER set for production), MC_OUT, REPO_ROOT.
suppressMessages({library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang)})

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
getc <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) x else d }
ROOT <- getc("REPO_ROOT", ".")
OUT  <- getc("MC_OUT", file.path(ROOT, "simulation", "output", "maxclust"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

owd <- setwd(ROOT); on.exit(setwd(owd), add = TRUE)   # steCI/R and scenarios.R by relative path
for (f in list.files(file.path(ROOT, "steCI", "R"), pattern = "[.]R$", full.names = TRUE)) source(f)
# scenarios.R gives make_scenarios(), DGM_FIXED and T_HOR; the two repo copies are byte-identical,
# so take whichever ships (the ARC zip carries simulation/scenarios.R). Fail loudly if neither is.
scen <- Filter(file.exists, file.path(ROOT, c(
  "simulation/scenarios.R",
  "simulation/scenarios.R")))
stopifnot(length(scen) >= 1L); source(scen[1])
setwd(owd)
stopifnot(is.function(dat_sim), is.function(ste_prepare), is.function(make_scenarios),
          is.list(DGM_FIXED))

# --- CELLS: the 27 n = 5000 scenarios, the ONLY sample size timereg bins (default max.clust=1000)
CELLS <- make_scenarios()
CELLS <- CELLS[CELLS$n == 5000L, , drop = FALSE]
CELLS <- CELLS[order(CELLS$scenario), , drop = FALSE]
rownames(CELLS) <- NULL
stopifnot(nrow(CELLS) == 27L)

REPS       <- geti("MC_REPS", 100L)     # replicates per cell (production)
CHUNK_REPS <- geti("MC_CHUNK", 10L)     # replicates per array task
REPS_OVR   <- geti("MC_REPS_OVERRIDE", NA_integer_)   # smoke only
R  <- if (is.na(REPS_OVR)) REPS else REPS_OVR
CR <- min(CHUNK_REPS, R)
SEED_BASE <- 700000L
stopifnot(R < 1000L)   # seeds are spaced 1000 apart per cell (sd_i below); guard cross-cell reuse

want <- as.integer(strsplit(getc("MC_CELLS", paste(seq_len(nrow(CELLS)), collapse = ":")), "[:,]")[[1]])
want <- want[!is.na(want)]
stopifnot(length(want) >= 1L, all(want >= 1L & want <= nrow(CELLS)))

tasks <- do.call(rbind, lapply(want, function(ci)
  data.frame(cell_idx = ci, chunk = seq_len(ceiling(R / CR)), stringsAsFactors = FALSE)))

tid <- geti("SLURM_ARRAY_TASK_ID", 1L)
cat(sprintf("[mc] task %d of %d (cells %s ; %d reps/cell, %d/chunk ; n=5000, max.clust default 1000)\n",
            tid, nrow(tasks), paste(want, collapse = ","), R, CR))
if (tid > nrow(tasks)) { cat("[mc] beyond grid -- nothing to do\n"); quit(save = "no") }
tk   <- tasks[tid, ]
cell <- CELLS[tk$cell_idx, ]

# Cell, reps and chunk all in the filename: a smoke run (R2-tag) cannot poison production, and a
# changed grid cannot silently overwrite a previous chunk.
outfile <- file.path(OUT, sprintf("mc_c%02d_%s_R%d_chunk%03d.rds",
                            tk$cell_idx, cell$scenario, R, tk$chunk))
if (file.exists(outfile)) { cat("[mc] already complete:", basename(outfile), "\n"); quit(save = "no") }
stopifnot(grepl("^mc_.*[.]rds$", basename(outfile)),
          identical(normalizePath(dirname(outfile), mustWork = TRUE),
                    normalizePath(OUT, mustWork = TRUE)))

idx <- ((tk$chunk - 1L) * CR) + seq_len(CR); idx <- idx[idx <= R]
stopifnot(length(idx) >= 1L)
cat(sprintf("[mc] cell %d (%s) alpha.0=%.2f gamma.L=%.1f | chunk %d, replicates %d..%d\n",
            tk$cell_idx, cell$scenario, cell$alpha.0, cell$gamma.L, tk$chunk, min(idx), max(idx)))

# The published cluster fit (ci_analytic.R lines 61-65), refit with and without the pin. NOT
# passing max.clust uses timereg's default (1000) = the published path; max.clust = NULL is the
# exact subject clustering. Only the robust variance changes; gamma is pin-invariant.
base_fml <- survival::Surv(time.new, time.stop.new, event) ~ const(A.baseline) + L.baseline
# Branch EXPLICITLY. NOT passing max.clust uses timereg's default (=1000, the published path);
# pin passes max.clust = NULL for the exact subject clustering. Do NOT build a do.call() arg list
# and set args$max.clust <- NULL: in R that DELETES the key, so aalen would run the default in
# BOTH arms and se_pinned == se_default by construction (a guaranteed-zero, meaningless run).
fit_cluster <- function(ds, w, pin) {
  if (pin)
    timereg::aalen(base_fml, data = ds, n.sim = 0, weights = w, robust = 1,
                   clusters = ds$id, max.clust = NULL)
  else
    timereg::aalen(base_fml, data = ds, n.sim = 0, weights = w, robust = 1,
                   clusters = ds$id)
}

gamma <- se_def <- se_pin <- rep(NA_real_, length(idx))
nclust <- rep(NA_integer_, length(idx))
fail  <- character(length(idx))
t0 <- Sys.time()
for (i in seq_along(idx)) {
  sd_i <- SEED_BASE + 1000L * tk$cell_idx + idx[i]   # collision-free per (cell, rep)
  r <- tryCatch({
    d  <- dat_sim(n = cell$n, n.visit = DGM_FIXED$n.visit, gamma.0 = cell$gamma.0,
                  gamma.L = cell$gamma.L, alpha.0 = cell$alpha.0, alpha.A = DGM_FIXED$alpha.A,
                  alpha.L = DGM_FIXED$alpha.L, alpha.U = DGM_FIXED$alpha.U, seed = sd_i)
    ds <- ste_prepare(d, n.visit = DGM_FIXED$n.visit)
    w  <- ds$ipw.stab.cum
    fd <- fit_cluster(ds, w, pin = FALSE)
    fp <- fit_cluster(ds, w, pin = TRUE)
    list(nc = length(unique(ds$id)), g = unname(fd$gamma[1, 1]), gp = unname(fp$gamma[1, 1]),
         sd = sqrt(fd$robvar.gamma[1, 1]), sp = sqrt(fp$robvar.gamma[1, 1]))
  }, error = function(e) { fail[i] <<- conditionMessage(e); NULL })
  if (is.null(r)) next
  # gamma MUST be pin-invariant (max.clust touches only the variance). Flag if a timereg change
  # ever breaks that, rather than silently averaging a corrupted point estimate.
  if (is.finite(r$g) && is.finite(r$gp) && abs(r$g - r$gp) > 1e-10 * max(1, abs(r$g)))
    fail[i] <- sprintf("gamma moved under pin: %.3e", r$g - r$gp)
  nclust[i] <- r$nc; gamma[i] <- r$g; se_def[i] <- r$sd; se_pin[i] <- r$sp
}
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

peak_rss_mb <- function() {
  if (!file.exists("/proc/self/status")) return(NA_real_)
  l <- grep("^VmHWM:", readLines("/proc/self/status", warn = FALSE), value = TRUE)
  if (!length(l)) return(NA_real_)
  as.numeric(sub("\\D*(\\d+).*", "\\1", l[1])) / 1024
}

# ATOMIC WRITE: temp in the SAME dir + rename, so a mid-write kill cannot leave a truncated .rds
# that the next submission skips as complete.
tmp <- paste0(outfile, ".part", Sys.getpid())
saveRDS(list(cell_idx = tk$cell_idx, scenario = cell$scenario, alpha.0 = cell$alpha.0,
             gamma.L = cell$gamma.L, gamma.0 = cell$gamma.0, n = cell$n, reps = R,
             chunk = tk$chunk, rep_idx = idx, seed = SEED_BASE + 1000L * tk$cell_idx + idx,
             nclust = nclust, gamma = gamma, se_default = se_def, se_pinned = se_pin, fail = fail,
             n_ok = sum(is.finite(se_def) & is.finite(se_pin)), secs = secs,
             rss_mb = peak_rss_mb()), tmp)
if (!file.rename(tmp, outfile)) { unlink(tmp); stop("could not rename ", tmp, " -> ", outfile) }

cat(sprintf("[mc] wrote %s | %d/%d ok | %.1fs (%.1fs/fit) | peak %s MB\n", basename(outfile),
            sum(is.finite(se_def)), length(idx), secs, secs / max(1L, 2L * length(idx)),
            ifelse(is.na(peak_rss_mb()), "n/a", sprintf("%.0f", peak_rss_mb()))))
if (any(nzchar(fail))) {
  cat(sprintf("[mc] WARNING: %d of %d replicates flagged\n", sum(nzchar(fail)), length(idx)))
  for (m in unique(fail[nzchar(fail)])) cat("   ", m, "\n")
}
