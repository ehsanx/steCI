# TRUNCATION EXPERIMENT -- ARC worker. One array task = one (cell x n x rep-chunk) -> one RDS.
#
# WHY. Our sharpest surviving novel claim -- "truncation restores calibration but destroys
# consistency: the untruncated bias SHRINKS with n while the truncated bias GROWS to 17% of the true
# effect at n=5000" -- rests on **R=50**, which `dev/probe_truncation_power.R` measured at **power
# 0.067** (z = 0.46). A coin flip. It must be redone at R that can actually see the effect (R=2000 =>
# power 0.83), and it must test whether the result is TRUNCATION or KEOGH'S CAP RULE.
#
# THE CAP RULE IS KEOGH'S, NOT OURS (analysis_seq.R:182; our port is bit-identical). It takes
# quantile(w, 0.95) over a distribution with a ~65% POINT MASS AT EXACTLY 1. Measured at
# gamma.L=3/gamma.0=-2: P(w==1) = 0.654, P(w<=1) = 0.959, so **the cap is EXACTLY 1.000** -- it deletes
# ALL up-weighting while keeping ALL down-weighting, and it FLIPS between 1.000 and >1 across
# replicates. That asymmetry is a candidate explanation for the toward-the-null bias.
#
# THE DESIGN EXPLOITS THE SHIPPED CODE: ci_lef(weight="trunc", truncate=q) passes q straight to
# ste_prepare, so the truncation QUANTILE is already a parameter. Varying q alone moves the cap ON and
# OFF the atom with no custom columns and no hand-rolled gradients:
#     stab        untruncated                     (what the simulation actually uses)
#     trunc@0.95  cap = 1.000 EXACTLY -> ON the atom      (Keogh's reference rule)
#     trunc@0.99  cap ~ 2.41          -> OFF the atom
#     trunc@0.999 cap ~ 10+           -> far off the atom
# If "grows with n" appears at 0.95 but NOT at 0.99/0.999, the result is the DISCONTINUITY.
# If it appears at all of them, it is truncation itself (and Crump 2009 predicts a bias, since
# trimming changes the estimand).
#
# ci_lef returns est/se/lower/upper, so this yields BIAS, CALIBRATION (se_ratio) and COVERAGE --
# the local laptop version could only do bias, which cannot test "restores calibration" at all.
#
#   Local smoke:  TRUNC_TASK=1 TRUNC_CHUNK=4 TRUNC_B=50 Rscript dev/truncation_arc.R
#   ARC:          see simulation/slurm/truncation.slurm
suppressWarnings(suppressMessages({library(dplyr); library(tidyr); library(timereg)
                                   library(survival); library(rlang); library(parallel)}))
geti <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) as.integer(x) else d }
getc <- function(v, d) { x <- Sys.getenv(v); if (nzchar(x)) x else d }

ROOT <- getc("REPO_ROOT", ".")
for (f in list.files(file.path(ROOT, "steCI/R"), pattern = "[.]R$", full.names = TRUE)) source(f)

TASK  <- geti("TRUNC_TASK", geti("SLURM_ARRAY_TASK_ID", 1L))
CHUNK <- geti("TRUNC_CHUNK", 32L)      # reps per task. MUST be <= cores: CHUNK=100 on 32 cores is the
                                        # bug that killed the main run at walltime (sim.slurm's own warning).
R_TOT <- geti("TRUNC_R", 2000L)
B     <- geti("TRUNC_B", 200L)
CORES <- geti("SLURM_CPUS_PER_TASK", max(1L, parallel::detectCores() - 1L))
OUT   <- getc("TRUNC_OUT", file.path(ROOT, "simulation/output/truncation"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
TAU <- 1:5

# The scientifically decisive trio, not a grid sweep:
#  strong_conf  -- the corner the R=50 claim was made in
#  gate_cell    -- gamma.L=0.5, gamma.0=0: BENIGN weights (ess_frac 0.941) yet a GENUINE ~12% variance
#                  shortfall (0.874 published, 0.884 Jensen-free). It has now defeated FIVE mechanisms
#                  AND the metric explanation. Truncation must be tested where the weights are fine.
#  clean_twin   -- gamma.L=0.5, gamma.0=-2: no shortfall at all (0.994). The negative control.
CELLS <- data.frame(
  label   = c("strong_conf", "gate_cell", "clean_twin"),
  gamma.L = c(3.0, 0.5, 0.5),
  gamma.0 = c(-2.0, 0.0, -2.0),
  alpha.0 = c(0.25, 0.25, 0.25), stringsAsFactors = FALSE)
NS <- c(300L, 1000L, 5000L)
REGIMES <- list(list(k = "stab",       weight = "stab",  truncate = 0.95),
                list(k = "trunc_q95",  weight = "trunc", truncate = 0.95),
                list(k = "trunc_q99",  weight = "trunc", truncate = 0.99),
                list(k = "trunc_q999", weight = "trunc", truncate = 0.999))

# Task -> (cell, n, chunk). n_chunks per (cell,n) = ceil(R_TOT / CHUNK).
NCH   <- ceiling(R_TOT / CHUNK)
COMBO <- expand.grid(chunk = seq_len(NCH), n = NS, ci = seq_len(nrow(CELLS)))
if (TASK > nrow(COMBO)) { cat("[trunc] task", TASK, "> ", nrow(COMBO), " -- nothing to do\n"); quit(save = "no") }
cmb <- COMBO[TASK, ]; cc <- CELLS[cmb$ci, ]; NN <- cmb$n
reps <- ((cmb$chunk - 1L) * CHUNK + 1L):min(cmb$chunk * CHUNK, R_TOT)
# R and CHUNK are IN THE FILENAME on purpose. Checkpointing skips an existing file, so without them a
# smoke run (R=8, chunk 001) and the full run (R=2000, chunk 001) would collide and the full run would
# silently REUSE the 8-rep smoke chunk. Distinct names make the smoke output inert instead of poisonous.
fout <- file.path(OUT, sprintf("trunc_%s_n%d_R%d_k%d_c%03d.rds", cc$label, NN, R_TOT, CHUNK, cmb$chunk))
if (file.exists(fout)) { cat("[trunc] exists, skipping:", fout, "\n"); quit(save = "no") }  # checkpointing

cat(sprintf("[trunc] task %d | %s | n=%d | chunk %d/%d | reps %d-%d | cores %d | B=%d\n",
            TASK, cc$label, NN, cmb$chunk, NCH, min(reps), max(reps), CORES, B))

one <- function(r) {
  seed <- 7e6 + cmb$ci * 1e5 + match(NN, NS) * 1e4 + r
  d <- dat_sim(n = NN, n.visit = 5, gamma.0 = cc$gamma.0, gamma.L = cc$gamma.L,
               alpha.0 = cc$alpha.0, alpha.A = -0.04, alpha.L = 0.015, alpha.U = 0.015, seed = seed)
  do.call(rbind, lapply(REGIMES, function(g) {
    lf <- try(ci_lef(d, t.hor = TAU, weight = g$weight, truncate = g$truncate,
                     B = B, approach = "a1", seed = seed + 1L), silent = TRUE)
    if (inherits(lf, "try-error")) return(NULL)
    m <- lf$mrd
    data.frame(regime = g$k, rep = r, n = NN, cell = cc$label,
               gamma.L = cc$gamma.L, gamma.0 = cc$gamma.0, alpha.0 = cc$alpha.0,
               tau = m$tau, method = m$type, estimate = m$estimate,
               se = m$se, lower = m$lower, upper = m$upper,
               stringsAsFactors = FALSE)
  }))
}

t0 <- proc.time()[["elapsed"]]
# NB the braces are load-bearing: at top level R closes `if (...) expr` at the newline and then chokes
# on a bare `else` ("unexpected 'else'"). mclapply forks (ARC/unix); Windows falls back to serial.
res <- if (.Platform$OS.type == "unix") {
  parallel::mclapply(reps, one, mc.cores = CORES)
} else {
  lapply(reps, one)
}
res <- do.call(rbind, res[!vapply(res, is.null, logical(1))])
if (is.null(res) || !nrow(res)) {
  cat("[trunc] WARNING: every replicate failed for task", TASK, "-- writing an empty marker\n")
  res <- data.frame()
}
saveRDS(res, fout)
cat(sprintf("[trunc] wrote %s : %d rows, %d reps, %.1f s\n",
            basename(fout), nrow(res),
            if (nrow(res)) length(unique(res$rep)) else 0L, proc.time()[["elapsed"]] - t0))
