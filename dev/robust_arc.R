# Re-run the three appendix robustness studies on ARC, SAVING PER-REPLICATE ESTIMATES.
#
# WHY THIS EXISTS. All three studies score raw coverage against ONE frozen Monte-Carlo draw of
# their own truth -- the same single-frozen-draw defect identified for the primary DGM. The truths are
# now available exactly (dev/robustness_truth_exact.R, dev/robustness_truth_cox.R), but the studies
# kept only 9 summary rows each and discarded the per-replicate estimates, so they cannot be
# re-scored. That, not the truth, is what forces a re-run.
#
# THE DESIGN POINT: THE TRUTH IS AN INPUT TO SCORING, NOT TO SIMULATION.
# This script computes NO truth and reads NO truth. It emits one row per replicate --
# (gamma_hat, se_model, se_row, se_cl) -- and nothing else. Scoring happens afterwards, on the
# workstation, in dev/robust_collect.R. Three consequences:
#   (1) tvfx's truth is still underived, and that does not block this run at all;
#   (2) any future truth correction is a seconds-long re-score, never another cluster job;
#   (3) the cluster never needs the exact-truth files, so they are not in the bundle.
#
# THE DGM AND FIT CODE IS NOT REIMPLEMENTED HERE. Each study's original script is read and
# source()d up to (but not including) its driver loop, so dat_sim2/dgm/fit2/fit_e1 are
# byte-identical to the code that produced the published tables. Reimplementing them would make
# any difference in the new numbers un-attributable: truth change, or transcription slip?
#
# SEEDING IS PRESERVED EXACTLY: seed = 1000*gi + r, with gi the ORIGINAL expand.grid row index and
# r the replicate index. So this reproduces the SAME replicates, and every change in a reported
# coverage is attributable to the truth alone. Do not "improve" the seeding.
#
# One array task = one (study, cell, replicate-chunk) -> one RDS. Checkpointed: an existing output
# is skipped, so a re-submit resumes in seconds.
#
# Env: ROBUST_STUDIES (colon-separated), ROBUST_REPS_OVERRIDE, ROBUST_OUT, REPO_ROOT
suppressMessages({library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang)})

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
getc <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) x else d }
ROOT <- getc("REPO_ROOT", ".")
OUT  <- getc("ROBUST_OUT", file.path(ROOT, "simulation", "output", "robust"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# --- study registry -------------------------------------------------------------------------
# Each study differs in three ways that MUST be dispatched, not assumed: the generator/fit
# function names, the grid's alpha column name, and the est matrix column names. Getting any of
# these wrong yields a silently mislabelled result rather than an error.
# Generator/fit are named by STRING, not captured as closures, because the study's definitions
# are loaded into a private environment (see below) and must be resolved from THERE. All three
# generators take a named `seed` argument, so one uniform call signature covers them.
STUDIES <- list(
  `2confounder` = list(
    file = "dev/robustness_2confounder.R", reps = 300L,
    ns = c(300L, 1000L), gen = "dat_sim2", fit = "fit2",
    cols = c("gamma", "model", "robust", "cluster")),
  misspec = list(
    file = "dev/robustness_misspec.R", reps = 500L,
    ns = c(300L, 1000L, 5000L), gen = "dgm", fit = "fit_e1",
    cols = c("gamma", "se_model", "se_row", "se_cl")),
  tvfx = list(
    file = "dev/robustness_tvfx.R", reps = 500L,
    ns = c(300L, 1000L, 5000L), gen = "dgm", fit = "fit_e1",
    cols = c("gamma", "se_model", "se_row", "se_cl")))

A0S <- c(0.15, 0.25, 0.40)

# Replicates per chunk, by n. Sized for ~15-25 min per ARC task from MEASURED laptop cost
# (0.84 / 2.50 / 30.4 s per replicate at n = 300 / 1000 / 5000) and the ~4x per-core slowdown
# measured on ARC during the truth run. Smaller chunks = better Slurm backfill; larger chunks =
# fewer files. Do NOT raise the n=5000 value without re-reading the walltime.
chunk_reps <- function(n) if (n >= 5000L) 10L else if (n >= 1000L) 100L else 250L

# The ORIGINAL grid, reproduced exactly. expand.grid varies the FIRST factor fastest, so
# gi = 1,2,3 is n = 300,1000,5000 at alpha.0 = 0.15. gi feeds the seed, so an ordering change
# here silently re-randomises every replicate while still "working".
study_grid <- function(s) expand.grid(n = STUDIES[[s]]$ns, alpha = A0S)

# --- flat task list -------------------------------------------------------------------------
want <- strsplit(getc("ROBUST_STUDIES", "2confounder:misspec:tvfx"), "[:,]")[[1]]
want <- want[nzchar(want)]
stopifnot(all(want %in% names(STUDIES)))

REPS_OVR <- geti("ROBUST_REPS_OVERRIDE", NA_integer_)   # smoke only; NEVER set for production
tasks <- do.call(rbind, lapply(want, function(s) {
  g <- study_grid(s); R <- if (is.na(REPS_OVR)) STUDIES[[s]]$reps else REPS_OVR
  do.call(rbind, lapply(seq_len(nrow(g)), function(gi) {
    cr <- min(chunk_reps(g$n[gi]), R)
    data.frame(study = s, gi = gi, n = g$n[gi], alpha = g$alpha[gi], reps = R,
               chunk = seq_len(ceiling(R / cr)), chunk_reps = cr, stringsAsFactors = FALSE)
  }))
}))

tid <- geti("SLURM_ARRAY_TASK_ID", 1L)
cat(sprintf("[robust] task %d of %d total tasks (studies: %s)\n",
            tid, nrow(tasks), paste(want, collapse = ",")))
if (tid > nrow(tasks)) { cat("[robust] beyond grid -- nothing to do\n"); quit(save = "no") }
tk <- tasks[tid, ]

# Study, cell, chunk and reps all in the filename: a smoke run must not be able to poison
# production, and a changed grid must not silently overwrite a previous run's chunk.
outfile <- file.path(OUT, sprintf("rob_%s_gi%02d_n%d_a%s_R%d_chunk%03d.rds",
                            tk$study, tk$gi, tk$n, sub("[.]", "", format(tk$alpha)),
                            tk$reps, tk$chunk))
if (file.exists(outfile)) { cat("[robust] already complete:", basename(outfile), "\n"); quit(save = "no") }

# --- load the study's ORIGINAL definitions, without its driver loop -------------------------
sf  <- file.path(ROOT, STUDIES[[tk$study]]$file)
src <- readLines(sf, warn = FALSE)
cut <- grep("^grid <- expand.grid", src)
stopifnot(length(cut) == 1L)      # if the driver line moved, FAIL rather than source the loop
tmp <- tempfile(fileext = ".R"); writeLines(src[1:(cut[1] - 1L)], tmp)

# SOURCE INTO A PRIVATE ENVIRONMENT, never the global one.
# The study scripts contain `for (f in list.files("steCI/R", ...)) source(f)`, whose loop
# variable `f` is exactly the name this script uses for its OUTPUT PATH. Sourcing globally let
# that loop clobber it, and saveRDS then wrote the results over steCI/R/weight_if.R -- a package
# source file. Caught by the local smoke run; git restored it. The lesson generalises beyond `f`:
# any name the study happens to use (src, cut, i, s, d, v, est, idx) could collide, so isolate
# the whole thing rather than renaming variables and hoping.
# parent = globalenv() so the study's functions still see the attached packages; their closures
# get ENV, so they resolve each other and their DGM constants normally.
ENV <- new.env(parent = globalenv())
owd <- setwd(ROOT); on.exit(setwd(owd), add = TRUE)   # the scripts source steCI/R by RELATIVE path
suppressMessages(sys.source(tmp, envir = ENV))
setwd(owd)

st  <- STUDIES[[tk$study]]
gen <- get(st$gen, envir = ENV)      # resolved from ENV, not from here
fit <- get(st$fit, envir = ENV)
stopifnot(is.function(gen), is.function(fit))

# Belt and braces. The private ENV should make this impossible, but the failure mode it guards
# against -- saveRDS() writing over a source file because a sourced script rebound the path
# variable -- is silent, destructive, and actually happened. Assert rather than trust.
stopifnot(is.character(outfile), length(outfile) == 1L,
          grepl("^rob_.*[.]rds$", basename(outfile)),
          identical(normalizePath(dirname(outfile), mustWork = TRUE),
                    normalizePath(OUT, mustWork = TRUE)))
idx <- ((tk$chunk - 1L) * tk$chunk_reps) + seq_len(tk$chunk_reps)
idx <- idx[idx <= tk$reps]
stopifnot(length(idx) >= 1L)

cat(sprintf("[robust] %s gi=%d n=%d alpha=%s | chunk %d, replicates %d..%d\n",
            tk$study, tk$gi, tk$n, format(tk$alpha), tk$chunk, min(idx), max(idx)))

est  <- matrix(NA_real_, length(idx), 4, dimnames = list(NULL, st$cols))
fail <- character(length(idx))
t0   <- Sys.time()
for (i in seq_along(idx)) {
  sd_i <- 1000L * tk$gi + idx[i]                   # the ORIGINAL seed formula -- do not change
  d <- tryCatch(gen(tk$n, tk$alpha, seed = sd_i),
                error = function(e) { fail[i] <<- paste("gen:", conditionMessage(e)); NULL })
  if (is.null(d)) next
  v <- tryCatch(fit(d), error = function(e) { fail[i] <<- paste("fit:", conditionMessage(e)); NULL })
  if (!is.null(v)) est[i, ] <- v
}
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# Peak RSS, so --mem is MEASURED on the cluster rather than estimated. Linux only.
peak_rss_mb <- function() {
  if (!file.exists("/proc/self/status")) return(NA_real_)
  l <- grep("^VmHWM:", readLines("/proc/self/status", warn = FALSE), value = TRUE)
  if (!length(l)) return(NA_real_)
  as.numeric(sub("\\D*(\\d+).*", "\\1", l[1])) / 1024
}

# ATOMIC WRITE. saveRDS() is not atomic, and this job is checkpointed on file EXISTENCE: a task
# killed by a walltime limit or a pre-emption mid-write would leave a truncated .rds that the
# next submission then skips as "already complete", silently baking a corrupt chunk into the
# results. Write to a temp name in the SAME directory (so the rename cannot cross a filesystem)
# and rename into place -- rename is atomic on POSIX, which is where this runs.
tmp_out <- paste0(outfile, ".part", Sys.getpid())
saveRDS(list(study = tk$study, gi = tk$gi, n = tk$n, alpha = tk$alpha, reps = tk$reps,
             chunk = tk$chunk, rep_idx = idx, seed = 1000L * tk$gi + idx,
             cols = st$cols, est = est, fail = fail,
             n_ok = sum(is.finite(est[, 1])), secs = secs, rss_mb = peak_rss_mb(),
             src_file = STUDIES[[tk$study]]$file, src_lines_sourced = cut[1] - 1L), tmp_out)
if (!file.rename(tmp_out, outfile)) {
  unlink(tmp_out)
  stop("could not rename ", tmp_out, " -> ", outfile)
}

cat(sprintf("[robust] wrote %s | %d/%d ok | %.1fs | peak %s MB\n", basename(outfile),
            sum(is.finite(est[, 1])), length(idx), secs,
            ifelse(is.na(peak_rss_mb()), "n/a", sprintf("%.0f", peak_rss_mb()))))
if (any(nzchar(fail))) {
  cat(sprintf("[robust] WARNING: %d of %d replicates failed\n", sum(nzchar(fail)), length(idx)))
  for (m in unique(fail[nzchar(fail)])) cat("   ", m, "\n")
}
