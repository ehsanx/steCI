# ADEMP simulation driver. One invocation = one SLURM array task = one
# (scenario x rep-chunk) cell. Runs locally too (defaults to task 1).
#
#   Local pilot:  Rscript simulation/run_sim.R           (task 1, tiny)
#   ARC array:    driven by simulation/slurm/sim.slurm via SLURM_ARRAY_TASK_ID
#
# Env vars: SLURM_ARRAY_TASK_ID, SLURM_CPUS_PER_TASK, SIM_N_REPS, SIM_CHUNK,
#           SIM_B, SIM_SEED_BASE, SIM_OUT, REPO_ROOT.

geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
getc <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) x else d }

REPO_ROOT <- getc("REPO_ROOT", ".")
suppressMessages({
  library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang)
})
# Load steCI: prefer the installed package, else source its R/ files found
# either under REPO_ROOT or as a sibling of the simulation/ folder.
load_steCI <- function(root) {
  if (requireNamespace("steCI", quietly = TRUE)) { library(steCI); return("installed") }
  cands <- unique(c(file.path(root, "steCI", "R"),
                    file.path(root, "..", "steCI", "R"),
                    "steCI/R", "../steCI/R"))
  hit <- cands[dir.exists(cands)][1]
  if (is.na(hit)) stop(
    "Cannot find the steCI package. Install it, or place the 'steCI/' folder ",
    "at REPO_ROOT (", normalizePath(root, mustWork = FALSE),
    ") or as a sibling of 'simulation/'. Searched: ", paste(cands, collapse = ", "),
    call. = FALSE)
  for (f in list.files(hit, pattern = "[.]R$", full.names = TRUE)) source(f)
  hit
}
load_steCI(REPO_ROOT)
source(file.path(REPO_ROOT, "simulation", "scenarios.R"))
source(file.path(REPO_ROOT, "simulation", "one_replicate.R"))

task    <- geti("SLURM_ARRAY_TASK_ID", 1L)
cores   <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
N_REPS  <- geti("SIM_N_REPS", 4L)      # per scenario (pilot default tiny)
CHUNK   <- geti("SIM_CHUNK", N_REPS)   # reps per array task
B       <- geti("SIM_B", 200L)         # bootstrap resamples
SEED_BASE <- geti("SIM_SEED_BASE", 20250228L)
OUT     <- getc("SIM_OUT", file.path(REPO_ROOT, "simulation", "output", "intermediate"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

scen <- make_scenarios()
n_chunks <- ceiling(N_REPS / CHUNK)
grid <- do.call(rbind, lapply(seq_len(nrow(scen)), function(i)
  data.frame(scen_id = scen$scen_id[i], chunk = seq_len(n_chunks))))
stopifnot(task >= 1L, task <= nrow(grid))

job <- grid[task, ]
cell <- scen[scen$scen_id == job$scen_id, ]
rep_from <- (job$chunk - 1L) * CHUNK + 1L
rep_to   <- min(job$chunk * CHUNK, N_REPS)
reps_idx <- rep_from:rep_to

# Shared replicate seeds (replicate i uses the same base seed in every scenario).
set.seed(SEED_BASE)
all_seeds <- sample.int(.Machine$integer.max, N_REPS)

JACK_MAXN <- geti("SIM_JACK_MAXN", 300L)  # set 0 to disable (jackknife is O(n) refits)
do_jack <- cell$n <= JACK_MAXN

cat(sprintf("[run_sim] task=%d scenario=%s (scen_id=%d) reps %d-%d B=%d cores=%d jack=%s\n",
            task, cell$scenario, cell$scen_id, rep_from, rep_to, B, cores, do_jack))

run_rep <- function(i) {
  # NOTE (2026-07-17): this tryCatch is now a BACKSTOP, not the failure mechanism. It used to be the
  # ONLY one, so any single method throwing (in practice .aalen_if/M6 going singular on a late risk
  # set with no A.baseline variation, which timereg::aalen survives) discarded the replicate for EVERY
  # method -- deletion conditioned on "M6 did not go singular", i.e. non-ignorable, and biased toward
  # our own recommended method. run_one() now isolates each method and records `failed_methods`, so
  # reaching this handler means the replicate failed OUTSIDE any method block (e.g. dat_sim itself).
  r <- tryCatch(
    run_one(cell, rep_seed = all_seeds[i], B = B, do_jackknife = do_jack),
    error = function(e) { message("rep ", i, " failed OUTSIDE any method block: ",
                                  conditionMessage(e)); NULL })
  if (is.null(r)) return(NULL)
  r$rep <- i; r$scenario <- cell$scenario; r$scen_id <- cell$scen_id
  r$n <- cell$n; r$alpha.0 <- cell$alpha.0
  r$gamma.0 <- cell$gamma.0; r$gamma.L <- cell$gamma.L
  r
}

outfile <- file.path(OUT, sprintf("sim_scen%03d_chunk%03d.rds", job$scen_id, job$chunk))
meta <- list(task = task, B = B, N_REPS = N_REPS, CHUNK = CHUNK, seed_base = SEED_BASE,
             R_version = as.character(getRversion()),
             timereg = as.character(utils::packageVersion("timereg")))

# --- Resume: reload any reps completed by a previous (possibly timed-out) run.
acc <- NULL; done <- integer(0)
if (file.exists(outfile)) {
  prev <- tryCatch(readRDS(outfile), error = function(e) NULL)
  # Refuse to resume across a settings change. Without this the block merges old replicates into
  # a new run whenever the file name matches, so a changed B, seed base or weight scheme would
  # produce a chunk that is half one configuration and half another -- undetectable downstream,
  # because nothing records per-replicate settings. The canonical logs contain no "resuming"
  # lines, so no shipped chunk is contaminated; this guards future users.
  if (!is.null(prev$meta)) {
    key <- function(m) c(B = m$B, N_REPS = m$N_REPS, CHUNK = m$CHUNK, seed_base = m$seed_base)
    if (!identical(key(prev$meta), key(list(B = B, N_REPS = N_REPS, CHUNK = CHUNK,
                                            seed_base = SEED_BASE))))
      stop("refusing to resume ", outfile, ": its meta does not match this run's settings.\n",
           "  Delete the file to start fresh, or restore the original settings.")
  }
  if (!is.null(prev$reps)) { acc <- prev$reps; done <- unique(prev$reps$rep) }
}
todo <- setdiff(reps_idx, done)
if (length(todo) == 0L) {
  cat(sprintf("[run_sim] task %d already complete (%d reps) -- skipping\n", task, length(done)))
  quit(save = "no", status = 0)
}
if (length(done)) cat(sprintf("[run_sim] resuming: %d done, %d to do\n", length(done), length(todo)))

is_unix <- .Platform$OS.type == "unix"
run_batch <- function(idx) {
  if (is_unix && cores > 1L) parallel::mclapply(idx, run_rep, mc.cores = cores)
  else lapply(idx, run_rep)
}

# --- Checkpoint after each batch of ~cores reps, so a walltime kill loses at
# most one batch. Write to a temp file then rename (never leave a half file).
batch_size <- max(1L, cores)
batches <- split(todo, ceiling(seq_along(todo) / batch_size))
for (bi in seq_along(batches)) {
  new <- do.call(rbind, Filter(Negate(is.null), run_batch(batches[[bi]])))
  acc <- if (is.null(acc)) new else rbind(acc, new)
  tmp <- paste0(outfile, ".tmp")
  saveRDS(list(reps = acc, cell = cell, meta = meta), tmp)
  file.rename(tmp, outfile)
  cat(sprintf("[run_sim] checkpoint %d/%d batches; %d/%d reps done -> %s\n",
              bi, length(batches), length(unique(acc$rep)), length(reps_idx), outfile))
}
cat("[run_sim] task", task, "complete:", nrow(acc), "rows\n")
