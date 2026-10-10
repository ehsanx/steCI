# M8 on ARC: re-estimate the E1 target `gamma_true` properly.
#
# WHY. The shipped target is ONE frozen Monte-Carlo draw. Measured replication SD is
# 0.00254 / 0.00382 / 0.00629 at alpha.0 = 0.15 / 0.25 / 0.4, against an E1 bias of ~0.0057 --
# i.e. the noise in the target is comparable to the quantity it is used to measure, and being a
# single draw shared across all replicates and scenarios it is a FIXED displacement that does
# not average out.
#
# TWO THINGS THIS FIXES, and the second is why the laptop plan was wrong in kind, not just in size:
#   (1) VARIANCE. Averaging K draws cuts the MC SE by sqrt(K). K=2000 -> ~0.00014, ~2.5% of the bias.
#   (2) n_hd BIAS. Averaging over seeds converges to E[gamma_hat(n_hd)], NOT to the population
#       gamma. If the aalen const() fit carries finite-sample bias at n_hd = 3e4, no amount of
#       seed-averaging removes it. So this script also sweeps n_hd, which is only affordable
#       because the job is on a cluster. A laptop budget silently forces (1) and ignores (2).
#
# One array task = one (alpha.0 x n_hd x seed-chunk) -> one RDS. Checkpointed: an existing output
# is skipped, so a re-submit resumes in seconds. Uses mclapply (fork) -- on Linux this needs no
# export, unlike the PSOCK path used on Windows.
#
# Env: TRUTH_K, TRUTH_CHUNK, TRUTH_NHD (comma-separated list), TRUTH_OUT, REPO_ROOT
suppressMessages({library(parallel); library(timereg); library(survival)})

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
getc <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) x else d }
ROOT <- getc("REPO_ROOT", ".")
for (f in list.files(file.path(ROOT, "steCI", "R"), pattern = "[.]R$", full.names = TRUE)) source(f)
source(file.path(ROOT, "simulation", "scenarios.R"))

K      <- geti("TRUTH_K", 2000L)
CHUNK  <- geti("TRUTH_CHUNK", 32L)
# Delimiter is a COLON, not a comma: sbatch --export is itself comma-separated, so a comma here
# silently truncates the list to its first element (this bit us on the first smoke run -- the
# sweep never ran and 6 of 9 tasks were no-ops). Commas still accepted for direct local use.
NHDS   <- as.integer(strsplit(getc("TRUTH_NHD", "30000"), "[:,]")[[1]])
NHDS   <- NHDS[is.finite(NHDS) & NHDS > 0]
stopifnot(length(NHDS) >= 1)
OUT    <- getc("TRUTH_OUT", file.path(ROOT, "simulation", "output", "truth"))
# TRUTH_CORES overrides the core count independently of the Slurm allocation. This is the lever for
# the large-n_hd arm: timereg::aalen is ~quadratic in n, so at n_hd = 3e5 each worker fits 6e5 rows,
# and 32 of those CONCURRENTLY is the plausible way this job dies. Lower it for that arm rather than
# raising --mem and hoping.
CORES  <- geti("TRUTH_CORES", geti("SLURM_CPUS_PER_TASK", 32L))

# Peak RSS, so the memory question is MEASURED on the cluster rather than estimated. Linux only;
# returns NA elsewhere. VmHWM is per-process high-water mark, so a forked child reports its own.
peak_rss_mb <- function() {
  f <- "/proc/self/status"
  if (!file.exists(f)) return(NA_real_)
  l <- grep("^VmHWM:", readLines(f, warn = FALSE), value = TRUE)
  if (!length(l)) return(NA_real_)
  as.numeric(sub("\\D*(\\d+).*", "\\1", l[1])) / 1024
}
# Slurm puts the task in a cgroup, whose peak covers ALL forked children together -- the number that
# actually decides whether --mem is enough.
# The first smoke run read 186 GB here -- that is the NODE's cgroup, not this job's, and is
# useless. Resolve the process's own cgroup path from /proc/self/cgroup first, and only fall back
# to the global files if that fails. Returns NA rather than a misleading number.
cgroup_peak_mb <- function() {
  paths <- character(0)
  if (file.exists("/proc/self/cgroup")) {
    ln <- readLines("/proc/self/cgroup", warn = FALSE)
    rel <- sub("^[0-9]+:[^:]*:", "", ln)
    rel <- rel[nzchar(rel) & rel != "/"]
    for (r in rel) paths <- c(paths,
      file.path("/sys/fs/cgroup", r, "memory.peak"),
      file.path("/sys/fs/cgroup/memory", r, "memory.max_usage_in_bytes"))
  }
  for (f in paths) {
    if (file.exists(f)) {
      v <- suppressWarnings(as.numeric(readLines(f, warn = FALSE)[1]))
      if (is.finite(v) && v > 0) return(v / 1024^2)
    }
  }
  NA_real_
}
A0S    <- c(0.15, 0.25, 0.4)
BASE   <- 91000000L
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

stopifnot(CHUNK <= CORES)   # CHUNK > cores is exactly the bug that cost the main run 10 cells

nchunk <- ceiling(K / CHUNK)
grid   <- expand.grid(chunk = seq_len(nchunk), nhd = NHDS, a0 = A0S)  # chunk varies fastest
tid    <- geti("SLURM_ARRAY_TASK_ID", 1L)
if (tid > nrow(grid)) { cat("[truth] task", tid, "> grid size", nrow(grid), "-- nothing to do\n"); quit(save = "no") }
g   <- grid[tid, ]
a0  <- g$a0; nhd <- g$nhd; ch <- g$chunk

# K and CHUNK and n_hd in the filename: a smoke run must not be able to poison production.
f <- file.path(OUT, sprintf("truth_a%s_nhd%d_K%d_C%d_chunk%03d.rds",
                            sub("[.]", "", format(a0)), nhd, K, CHUNK, ch))
if (file.exists(f)) { cat("[truth] already complete:", basename(f), "\n"); quit(save = "no") }

# Globally unique seeds: no (a0, nhd, chunk, draw) combination can collide.
lvl  <- match(a0, A0S); hix <- match(nhd, NHDS)
idx  <- ((ch - 1L) * CHUNK) + seq_len(CHUNK)
idx  <- idx[idx <= K]
seed <- BASE + (lvl - 1L) * 10000000L + (hix - 1L) * 1000000L + idx

cat(sprintf("[truth] task %d: alpha.0=%s n_hd=%d chunk %d/%d draws %d..%d on %d cores\n",
            tid, format(a0), nhd, ch, nchunk, min(idx), max(idx), CORES))
t0 <- Sys.time()
res <- mclapply(seq_along(idx), function(i) {
  tv <- try(true_values(n = nhd, n.visit = DGM_FIXED$n.visit, t.hor = T_HOR,
                        alpha.0 = a0, alpha.A = DGM_FIXED$alpha.A,
                        alpha.L = DGM_FIXED$alpha.L, alpha.U = DGM_FIXED$alpha.U,
                        seed = seed[i], n_hd = nhd), silent = TRUE)
  list(gamma = if (inherits(tv, "try-error")) NA_real_ else tv$gamma_true,
       rss   = peak_rss_mb())
}, mc.cores = CORES)

# mclapply returns a try-error OBJECT (not an element of our list) if a child dies outright, e.g. is
# OOM-killed. Handle that, or an OOM would surface as an unhelpful subscript error.
ok      <- vapply(res, function(z) is.list(z) && !is.null(z$gamma), logical(1))
g_vals  <- rep(NA_real_, length(res)); g_vals[ok] <- vapply(res[ok], function(z) z$gamma, numeric(1))
rss_all <- rep(NA_real_, length(res)); rss_all[ok] <- vapply(res[ok], function(z) z$rss, numeric(1))
n_died  <- sum(!ok)
secs    <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
rss_max <- suppressWarnings(max(rss_all, na.rm = TRUE)); if (!is.finite(rss_max)) rss_max <- NA_real_
cg      <- cgroup_peak_mb()

saveRDS(list(alpha0 = a0, n_hd = nhd, chunk = ch, K = K, CHUNK = CHUNK, cores = CORES,
             draw_idx = idx, seed = seed, gamma = g_vals,
             n_ok = sum(is.finite(g_vals)), n_died = n_died,
             rss_per_worker_mb = rss_all, rss_max_mb = rss_max, cgroup_peak_mb = cg,
             secs = secs), f)

cat(sprintf("[truth] wrote %s | %d/%d finite | mean %.6f | %.1fs\n",
            basename(f), sum(is.finite(g_vals)), length(idx), mean(g_vals, na.rm = TRUE), secs))
if (n_died) cat(sprintf("[truth] WARNING: %d of %d workers died outright (OOM kill looks like this)\n",
                        n_died, length(res)))
# The line to read after the smoke run: it decides whether the n_hd = 3e5 arm fits in --mem.
cat(sprintf("[truth] MEMORY n_hd=%d cores=%d | peak/worker %.0f MB | task cgroup peak %s MB | projected at 32 workers %s MB\n",
            nhd, CORES, rss_max,
            ifelse(is.na(cg), "n/a", sprintf("%.0f", cg)),
            ifelse(is.na(rss_max), "n/a", sprintf("%.0f", rss_max * 32))))
