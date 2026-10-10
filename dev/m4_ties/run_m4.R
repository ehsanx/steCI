# Runner for the M4 tie experiment (pilot and production). Resumable; one .rds per task.
#
#   Rscript dev/m4_ties/run_m4.R <stage>          stage in {pilot, prod, prod5000, rerun...}
#   env: M4_WORKERS (default 6)
#
# The task list for each stage is FIXED in tasks_<stage>.rds (written by make_tasks.R before the
# stage starts), so a restart resumes the same design: tasks whose output file exists are skipped.
# Workers source the code tree SNAP (env M4_SNAPSHOT; default the repository root). The original run
# used a frozen copy of the package code, so parallel work could not race it.
# PSOCK cluster (no fork on Windows).
suppressMessages(library(parallel))
ROOT <- normalizePath(Sys.getenv("REPO_ROOT", "."))
HERE <- file.path(ROOT, "dev/m4_ties")
SNAP <- normalizePath(Sys.getenv("M4_SNAPSHOT", ROOT))
args <- commandArgs(TRUE)
stage <- if (length(args)) args[1] else stop("usage: run_m4.R <stage>")
NW <- as.integer(Sys.getenv("M4_WORKERS", "6"))
tasks <- readRDS(file.path(HERE, sprintf("tasks_%s.rds", stage)))
OUT <- file.path(HERE, "reps", stage)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
tasks$file <- file.path(OUT, sprintf("s%03d_r%04d_%s.rds", tasks$scen_id, tasks$rep, tasks$mode))
done <- file.exists(tasks$file)
# M4_REUSE=stage1,stage2: also skip a task whose output already exists under another stage (the
# rerun reuses the production replicates, which are the same canonical replicates).
for (rs in strsplit(Sys.getenv("M4_REUSE"), ",")[[1]]) if (nzchar(rs))
  done <- done | file.exists(file.path(HERE, "reps", rs, basename(tasks$file)))
todo <- tasks[!done, ]
# M4_MAXREP=24: restrict this runner to replicates <= 24 (a helper that brings the interim forward).
if (nzchar(Sys.getenv("M4_MAXREP"))) todo <- todo[todo$rep <= as.integer(Sys.getenv("M4_MAXREP")), ]
# M4_ONLY_N=5000: restrict this runner to one sample size (a helper runner for the slowest part).
if (nzchar(Sys.getenv("M4_ONLY_N")))
  todo <- todo[c(300L, 1000L, 5000L)[((todo$scen_id - 1L) %% 3L) + 1L] == as.integer(Sys.getenv("M4_ONLY_N")), ]
# M4_REVERSE=1: work the list from the end, so a second runner started on the same stage (to add
# workers) meets the first in the middle; a task both pick up is computed twice, identically
# (deterministic), and the atomic write makes the duplicate harmless.
if (identical(Sys.getenv("M4_REVERSE"), "1")) todo <- todo[rev(seq_len(nrow(todo))), ]
LOG <- file.path(HERE, sprintf("run_%s.log", stage))
cat(sprintf("[%s] stage=%s tasks=%d todo=%d workers=%d\n", format(Sys.time()), stage,
            nrow(tasks), nrow(todo), NW), file = LOG, append = TRUE)
if (!nrow(todo)) quit(save = "no")

cl <- makePSOCKcluster(NW)
on.exit(stopCluster(cl), add = TRUE)
clusterExport(cl, c("HERE", "SNAP"))
clusterEvalQ(cl, {
  suppressMessages({ library(dplyr); library(tidyr); library(timereg); library(survival); library(rlang) })
  for (f in list.files(file.path(SNAP, "steCI/R"), "[.]R$", full.names = TRUE)) source(f)
  source(file.path(SNAP, "simulation/scenarios.R"))
  source(file.path(HERE, "tie_exact.R")); source(file.path(HERE, "prepare_fast.R"))
  source(file.path(HERE, "m4_core.R"))
  SCEN <- make_scenarios()
  set.seed(20250228L); ALL_SEEDS <- sample.int(.Machine$integer.max, 1000L)   # run_sim.R, N_REPS=1000
  NFAIL <- readRDS(file.path(HERE, "canon_nfail.rds"))
  NULL
})
one_task <- function(k, todo, LOG) {
  tk <- todo[k, ]
  if (file.exists(tk$file)) return("skip")
  cell <- SCEN[SCEN$scen_id == tk$scen_id, ]
  nf <- NFAIL$boot_n_fail[NFAIL$scen_id == tk$scen_id & NFAIL$rep == tk$rep]
  nf <- if (length(nf)) nf else NA_integer_
  res <- tryCatch(paired_rep(cell, tk$rep, ALL_SEEDS[tk$rep], mode = tk$mode, B = 200L,
                             canon_nfail = nf),
                  error = function(e) list(status = paste("ERROR:", conditionMessage(e)),
                                           cell = cell$scenario, rep = tk$rep))
  # an advance-mode replicate whose prepare failures do not match the canonical boot_n_fail had a
  # canonical refit fail INSIDE aalen(), so the RNG advance cannot be trusted: redo it in full mode
  # (README.md here). Saved under the task's file name with mode = "full" recorded inside.
  if (identical(res$status, "ok") && identical(res$mode, "advance") && isFALSE(res$advance_aligned))
    res <- tryCatch(paired_rep(cell, tk$rep, ALL_SEEDS[tk$rep], mode = "full", B = 200L, canon_nfail = nf),
                    error = function(e) list(status = paste("ERROR:", conditionMessage(e)),
                                             cell = cell$scenario, rep = tk$rep))
  tmp <- paste0(tk$file, ".tmp")
  saveRDS(res, tmp); file.rename(tmp, tk$file)
  cat(sprintf("[%s] s%03d r%04d %-7s %-30s %s %.0fs\n", format(Sys.time(), "%H:%M:%S"), tk$scen_id,
              tk$rep, tk$mode, cell$scenario, res$status,
              if (is.null(res$seconds)) NA else res$seconds), file = LOG, append = TRUE)
  res$status
}
st <- parLapplyLB(cl, seq_len(nrow(todo)), one_task, todo = todo, LOG = LOG, chunk.size = 1L)
cat(sprintf("[%s] stage=%s finished: %s\n", format(Sys.time()), stage,
            paste(names(table(unlist(st))), table(unlist(st)), collapse = ", ")), file = LOG, append = TRUE)
