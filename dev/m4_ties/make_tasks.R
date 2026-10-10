# Fixed task lists for each stage of the M4 tie experiment (written once, before the stage runs).
#   Rscript dev/m4_ties/make_tasks.R pilot
# The production task lists are written from the fixed design (see README.md here and the prod branch below).
ROOT <- Sys.getenv("REPO_ROOT", ".")
SNAP <- Sys.getenv("M4_SNAPSHOT", ROOT)
source(file.path(SNAP, "simulation/scenarios.R"))
HERE <- file.path(ROOT, "dev/m4_ties")
s <- make_scenarios()
stage <- commandArgs(TRUE)[1]
id <- function(n, g0 = NULL, gl = NULL, a0 = NULL) {
  k <- s$n == n
  if (!is.null(g0)) k <- k & s$gamma.0 %in% g0
  if (!is.null(gl)) k <- k & s$gamma.L %in% gl
  if (!is.null(a0)) k <- k & s$alpha.0 %in% a0
  s$scen_id[k]
}
tk <- function(ids, reps, mode) expand.grid(scen_id = ids, rep = reps, mode = mode,
                                            stringsAsFactors = FALSE)
if (stage == "pilot") {
  # PILOT: canonical replicates 999 and 1000 only (never used by production, which uses 1..K).
  t <- rbind(
    tk(id(5000L, -1, 1.5), 1000L, "advance"),              # n=5000 timing, 3 cells
    tk(id(1000L, -1, 1.5), 1000L, "full"),                 # n=1000 full-mode validation, 3 cells
    tk(id(300L), c(999L, 1000L), "full"),                  # n=300 all 27 cells, full mode
    tk(id(1000L), c(999L, 1000L), "advance"),              # n=1000 all 27 cells, advance mode
    data.frame(scen_id = c(4L, 7L, 7L, 16L, 16L, 25L),     # canonical reps with failed resamples:
               rep = c(1000L, 999L, 1000L, 999L, 1000L, 999L),   # advance-mode alignment check
               mode = "advance", stringsAsFactors = FALSE))
} else if (stage == "prod") {
  # PRODUCTION: canonical replicates 1..K, disjoint from the pilot.
  K <- as.integer(commandArgs(TRUE)[2:4]); names(K) <- c("300", "1000", "5000")
  stopifnot(!anyNA(K), all(K >= 24L))
  NFULL <- c("300" = 1L, "1000" = 1L)                   # full-mode validation subset
  c1000_full <- id(1000L, g0 = -1, gl = 1.5)            # 3 cells (medium prevalence, moderate confounding)
  c5000 <- id(5000L, gl = 1.5)                          # 9 cells: 3 event rates x 3 prevalences
  t <- rbind(
    tk(c5000, seq_len(K["5000"]), "advance"),
    tk(c1000_full, seq_len(NFULL["1000"]), "full"),
    tk(id(300L), seq_len(NFULL["300"]), "full"),
    tk(id(1000L), seq_len(K["1000"]), "advance"),
    tk(id(300L), seq_len(K["300"]), "advance"))
  # a replicate run in full mode is not also run in advance mode
  fullkey <- paste(t$scen_id, t$rep)[t$mode == "full"]
  t <- t[!(t$mode == "advance" & paste(t$scen_id, t$rep) %in% fullkey), ]
  # interleave so every cell accrues replicates at the same pace (partial results stay balanced)
  t <- t[order(t$rep, t$mode != "full", s$n[match(t$scen_id, s$scen_id)], t$scen_id), ]
} else if (stage == "rerun") {
  # FULL tie-exact M4 rerun on the canonical resamples: all 81 cells x canonical replicates 1..1000,
  # advance mode (E1; E2 on flagged resamples). Run with M4_REUSE=prod so production replicates are
  # not recomputed; assemble with  Rscript dev/m4_ties/assemble_rerun.R prod rerun
  t <- tk(s$scen_id, 1:1000, "advance")
  t <- t[order(t$rep, s$n[match(t$scen_id, s$scen_id)], t$scen_id), ]
} else stop("unknown stage ", stage)
t <- t[!duplicated(t), ]
saveRDS(t, file.path(HERE, sprintf("tasks_%s.rds", stage)))
cat(stage, ":", nrow(t), "tasks\n"); print(table(s$n[match(t$scen_id, s$scen_id)], t$mode))
