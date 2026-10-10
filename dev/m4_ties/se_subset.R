# Truth-free descriptive summary of the COMPLETED M4 tie-experiment replicates (reduced scope).
#
#   Rscript dev/m4_ties/se_subset.R   -> se_subset.rds, se_subset.csv, se_subset_table.tex + printout
#
# Scope: the production run was stopped at 717 of 14,526 tasks and was not completed. The
# pre-specified decision rule is not evaluated and no coverage is computed. This script never opens dev/truth_cache_exact.rds: it
# compares the canonical jitter arm (timereg, canon_boot.rds) with the tie-exact arm (E2, else E1)
# on the SAME canonical replicate and the SAME B = 200 resamples, for the bootstrap SE and the
# percentile interval only. The last block reads the PUBLISHED, truth-free per-scenario SE ratios
# (canonical metrics.rds, column se_ratio only) to express the SE change in Table e1's units.
#
# Units: a "pair" is one (scenario, replicate). The pilot (replicates 999-1000) and the production
# subset (replicates 1-24) are disjoint and are summarised SEPARATELY, never pooled. Full-mode
# validation re-runs that duplicate an advance-mode pair are dropped (one copy kept; the exact arms
# agree to 1e-14). Primary production subset = replicate indices complete in every scenario of that
# n (balanced); the incomplete indices (n=1000 rep 11, n=5000 rep 20) appear only in the
# "all completed" sensitivity rows.
#
# Summaries per (n, estimand): SE change = mean over scenarios of the per-scenario ratio of mean SEs
# (timereg / tie-exact) - 1; percentile-width change likewise; lower/upper percentile-limit shifts
# = mean of (timereg - tie-exact) / tie-exact SE. Monte-Carlo SEs: replicate i uses the same seed in
# every scenario (run_sim.R), so pairs sharing a replicate index are dependent; MCSEs are a
# delete-one-replicate-index jackknife (K = 11/10/14 indices). With K <= 2 (pilot) no MCSE is given.
suppressMessages(library(stats))
#
# Inputs: canon_boot.rds (32 MB; rebuild with extract_canonical.R from data/exact/results/
# replicates_all.rds) and reps/{pilot,prod}/ (837 replicate files, 15 MB) are NOT in git (see the
# README, "Large files"); M4_DATA names the directory holding them (default dev/m4_ties). Outputs
# go to M4_OUT (default dev/m4_ties): se_subset.rds/.csv/_table.tex as shipped here.
ROOT <- Sys.getenv("REPO_ROOT", ".")
HERE <- file.path(ROOT, "dev/m4_ties")
DATA <- Sys.getenv("M4_DATA", HERE)
OUTD <- Sys.getenv("M4_OUT", HERE)
source(file.path(ROOT, "steCI/R/utils.R"))          # .boot_ci (the canonical summariser)
source(file.path(ROOT, "simulation/scenarios.R"))
SCEN <- make_scenarios()
EST <- c("gamma", paste0("mrd", 1:5))

cb <- readRDS(file.path(DATA, "canon_boot.rds"))
stopifnot(!any(grepl("truth|cov", names(cb), ignore.case = TRUE)))
cb <- cb[cb$method == "boot_percentile", ]
cb$key <- paste(cb$scen_id, cb$rep, cb$estimand, ifelse(is.na(cb$tau), "NA", cb$tau))

load_stage <- function(stage) {
  f <- list.files(file.path(DATA, "reps", stage), "[.]rds$", full.names = TRUE)
  R <- lapply(f, readRDS)
  R <- R[vapply(R, `[[`, "", "status") == "ok"]
  key <- paste(vapply(R, `[[`, 0L, "scen_id"), vapply(R, `[[`, 0L, "rep"))
  mode <- vapply(R, `[[`, "", "mode")
  o <- order(key, mode != "advance")                       # keep the advance copy of a duplicate
  R <- R[o][!duplicated(key[o])]
  rows <- list()
  for (r in R) {
    cell <- SCEN[SCEN$scen_id == r$scen_id, ]
    for (k in seq_along(EST)) {
      est <- if (k == 1) "hd" else "mrd"; tau <- if (k == 1) NA else k - 1
      cj <- cb[match(paste(r$scen_id, r$rep, est, ifelse(is.na(tau), "NA", tau)), cb$key), ]
      if (is.na(cj$se)) next
      xe <- ifelse(!is.na(r$E2[, k]), r$E2[, k], r$E1[, k])
      be <- .boot_ci(xe, cj$estimate, 0.95)[1, ]           # percentile row; se is type-free
      rows[[length(rows) + 1]] <- data.frame(
        stage = stage, scen_id = r$scen_id, n = cell$n, rep = r$rep, mode = r$mode,
        estimand = if (k == 1) "gamma" else "mrd", tau = tau,
        se_j = cj$se, se_e = be$se, lo_j = cj$lower, up_j = cj$upper, lo_e = be$lower, up_e = be$upper,
        nval_e = be$n_valid, nties_mean = mean(r$info[, "nties"], na.rm = TRUE),
        # out-of-range flag (MRD only): a tie-exact resample estimate, or a timereg percentile limit,
        # outside [-1, 1]. Reported as a sensitivity row; never dropped from the primary summaries.
        oor = k > 1 && (any(abs(xe) > 1, na.rm = TRUE) || cj$lower < -1 || cj$upper > 1),
        stringsAsFactors = FALSE)
    }
  }
  d <- do.call(rbind, rows)
  validse <- function(se, est) is.finite(se) & se > 0 & (est != "mrd" | se <= 0.5)
  d$ok <- validse(d$se_j, d$estimand) & validse(d$se_e, d$estimand)
  d
}

P <- rbind(load_stage("pilot"), load_stage("prod"))
P$w_j <- P$up_j - P$lo_j; P$w_e <- P$up_e - P$lo_e
P$infl <- P$se_j / P$se_e - 1                                   # per-pair SE change (timereg/exact - 1)
P$winfl <- P$w_j / P$w_e - 1                                    # per-pair percentile-width change
P$dlo <- (P$lo_j - P$lo_e) / P$se_e                             # lower-limit shift, tie-exact SE units
P$dup <- (P$up_j - P$up_e) / P$se_e                             # upper-limit shift
P$dmid <- (P$dlo + P$dup) / 2                                   # midpoint shift

# balanced production subset: replicate indices present in every scenario of that n
P$balanced <- TRUE
for (nn in unique(P$n[P$stage == "prod"])) {
  s <- P$stage == "prod" & P$n == nn
  ncell <- length(unique(P$scen_id[s]))
  cpr <- tapply(P$scen_id[s], P$rep[s], function(x) length(unique(x)))
  P$balanced[s] <- P$rep[s] %in% as.integer(names(cpr)[cpr == ncell])
}

# ---- statistics and their replicate-index jackknife ------------------------------------------
rom <- function(d, a, b) mean(vapply(split(d, d$scen_id), function(x) mean(x[[a]]) / mean(x[[b]]), 0)) - 1
cmean <- function(d, col) mean(tapply(d[[col]], d$scen_id, mean))
jk <- function(d, f) {
  reps <- sort(unique(d$rep)); K <- length(reps)
  if (K <= 2) return(NA_real_)
  th <- vapply(reps, function(r) f(d[d$rep != r, ]), 0)
  sqrt((K - 1) / K * sum((th - mean(th))^2))
}
summ <- function(d) {
  d <- d[d$ok, ]
  fse <- function(x) rom(x, "se_j", "se_e"); fw <- function(x) rom(x, "w_j", "w_e")
  flo <- function(x) cmean(x, "dlo"); fup <- function(x) cmean(x, "dup"); fmid <- function(x) cmean(x, "dmid")
  dd <- d[!d$oor, ]
  data.frame(cells = length(unique(d$scen_id)), K = length(unique(d$rep)), pairs = nrow(d),
             se_infl = fse(d), se_infl_mcse = jk(d, fse),
             se_infl_median = median(d$infl), se_infl_max = max(d$infl),
             share_pos = mean(d$infl > 0), share_gt5 = mean(d$infl > 0.05),
             w_infl = fw(d), w_infl_mcse = jk(d, fw),
             dlo = flo(d), dlo_mcse = jk(d, flo), dup = fup(d), dup_mcse = jk(d, fup),
             dmid = fmid(d), dmid_mcse = jk(d, fmid),
             n_oor = sum(d$oor), se_infl_excl_oor = if (nrow(dd)) fse(dd) else NA)
}
grp <- function(D, label) {
  out <- do.call(rbind, lapply(split(D, list(D$n, D$estimand, ifelse(is.na(D$tau), 0, D$tau)), drop = TRUE),
    function(d) cbind(data.frame(set = label, n = d$n[1], estimand = d$estimand[1], tau = d$tau[1]), summ(d))))
  out[order(out$n, out$estimand != "gamma", out$tau), ]
}
S <- rbind(grp(P[P$stage == "pilot", ], "pilot"),
           grp(P[P$stage == "prod" & P$balanced, ], "production (balanced)"),
           grp(P[P$stage == "prod", ], "production (all completed)"))
rownames(S) <- NULL

# ---- implied Table e1 M4 SE-ratio column (truth-free) ------------------------------------------
# The published SE ratio (mean SE / empirical SD) is truth-free, and the empirical SD does not
# depend on the tie convention (the original datasets have no ties). Rescaling each scenario's
# published M4 (percentile) SE ratio by that scenario's paired ratio of means mean(se_e)/mean(se_j)
# estimates the value tie-exact refits would give -- the construction of decide_m4.R's se_shift_n,
# applied to the balanced production subset (all 27 scenarios at n = 300 and 1000, from 11 and 10 of
# the 1000 replicates). n = 5000 covers 9 of the row's 27 scenarios and is shown for that subset only.
met <- readRDS(file.path(ROOT, "data/exact/results/metrics.rds"))
cimet <- met$ci[met$ci$method == "boot_percentile" & met$ci$estimand == "hd", c("scenario", "se_ratio")]
rm(met)
cimet$scen_id <- SCEN$scen_id[match(cimet$scenario, SCEN$scenario)]
te1 <- readRDS(file.path(HERE, "canonical_rows.rds"))$tab_e1
te1 <- te1[te1$method == "boot_percentile", ]
E1 <- do.call(rbind, lapply(c(300L, 1000L, 5000L), function(nn) {
  d <- P[P$stage == "prod" & P$balanced & P$n == nn & P$estimand == "gamma" & P$ok, ]
  d$c <- cimet$se_ratio[match(d$scen_id, cimet$scen_id)]; stopifnot(!anyNA(d$c))
  f <- function(x) mean(vapply(split(x, x$scen_id), function(g) g$c[1] * (mean(g$se_e) / mean(g$se_j) - 1), 0))
  v0 <- te1$se_ratio[te1$n == nn]; sh <- f(d); ns <- length(unique(d$scen_id))
  data.frame(n = nn, scenarios = ns, K = length(unique(d$rep)), published_row = v0,
             shift = sh, mcse = jk(d, f), implied_row = if (ns == 27) v0 + sh else NA,
             printed_now = sprintf("%.2f", v0),
             printed_tie_exact = if (ns == 27) sprintf("%.2f", v0 + sh) else "n/a (9 of 27 scenarios)")
}))

saveRDS(list(pairs = P, summary = S, e1_ser = E1, note = "truth-free; no coverage; see header"),
        file.path(OUTD, "se_subset.rds"))
write.csv(S, file.path(OUTD, "se_subset.csv"), row.names = FALSE)

# ---- printout ---------------------------------------------------------------------------------
op <- options(width = 220, digits = 3)
up <- function(s) nrow(unique(P[s, c("scen_id", "rep")]))
cat(sprintf("unique pairs: pilot %d ; prod %d ; prod balanced %d\n", up(P$stage == "pilot"),
            up(P$stage == "prod"), up(P$stage == "prod" & P$balanced)))
cat(sprintf("scenarios: pilot %d ; prod %d ; invalid-SE rows dropped: %d\n",
            length(unique(P$scen_id[P$stage == "pilot"])), length(unique(P$scen_id[P$stage == "prod"])), sum(!P$ok)))
for (st in c("pilot", "prod")) {
  x <- P[P$stage == st & P$ok & (st == "pilot" | P$balanced), ]
  cat(sprintf("%s: timereg SE larger in %d of %d pair-estimands (%.2f%%); per-row minimum share %.3f\n", st,
              sum(x$infl > 0), nrow(x), 100 * mean(x$infl > 0),
              min(S$share_pos[S$set == if (st == "pilot") "pilot" else "production (balanced)"])))
}
cat(sprintf("balanced K by n: %s ; pilot n=5000 scenarios: %s\n",
            paste(tapply(P$rep[P$stage == "prod" & P$balanced], P$n[P$stage == "prod" & P$balanced],
                         function(x) length(unique(x))), collapse = "/"),
            paste(sort(unique(P$scen_id[P$stage == "pilot" & P$n == 5000])), collapse = " ")))
cat("\n"); print(S, row.names = FALSE)
cat("\n== implied Table e1 M4 SE ratio under tie-exact refits (balanced production subset) ==\n")
print(E1, row.names = FALSE, digits = 4)
options(op)

# ---- LaTeX (Web Table tab:boottie) ---------------------------------------------------------------
# Panel A (production, balanced): n | estimand | SE change % (MCSE) | median % | % of replicates > 5%
#   | width change % (MCSE) | lower-limit shift (MCSE) | upper-limit shift (MCSE)
# Panel B (pilot): estimand x n, mean SE change % (no MCSE; K <= 2).
pc <- function(x) ifelse(is.na(x), "--", ifelse(x > 0 & 100 * x < 0.05, "$<$0.1", sprintf("%.1f", 100 * x)))
pm <- function(x) ifelse(is.na(x), "--", ifelse(100 * x < 0.05, "$<$0.1", sprintf("%.1f", 100 * x)))
sg <- function(x) ifelse(is.na(x), "--", sub("^-", "$-$", sprintf("%.2f", ifelse(round(x, 2) == 0, 0, x))))
sm <- function(x) ifelse(is.na(x), "--", ifelse(x < 0.005, "$<$0.01", sprintf("%.2f", x)))
lab <- function(e, t) ifelse(e == "gamma", "$\\gamma$", sprintf("$\\tau=%d$", t))
Pb <- S[S$set == "production (balanced)", ]; Pp <- S[S$set == "pilot", ]
A <- vapply(seq_len(nrow(Pb)), function(i) { a <- Pb[i, ]
  sprintf("%s & %s & %s (%s) & %s & %.0f & %s (%s) & %s (%s) & %s (%s) \\\\",
    if (i == 1 || Pb$n[i - 1] != a$n) sprintf("$%d$", a$n) else "", lab(a$estimand, a$tau),
    pc(a$se_infl), pm(a$se_infl_mcse), pc(a$se_infl_median), 100 * a$share_gt5,
    pc(a$w_infl), pm(a$w_infl_mcse), sg(a$dlo), sm(a$dlo_mcse), sg(a$dup), sm(a$dup_mcse)) }, "")
B <- vapply(split(Pp, factor(paste(Pp$estimand, Pp$tau), levels = unique(paste(Pp$estimand, Pp$tau)))),
  function(x) sprintf(" & %s & %s & %s & %s & & & \\\\", lab(x$estimand[1], x$tau[1]),
                      pc(x$se_infl[x$n == 300]), pc(x$se_infl[x$n == 1000]), pc(x$se_infl[x$n == 5000])), "")
writeLines(c("% Panel A rows (production, balanced)", A, "% Panel B rows (pilot)", B),
           file.path(OUTD, "se_subset_table.tex"))
