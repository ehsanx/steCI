# Approach-2 vs Approach-1: % change in the mean MRD standard error, on COMMON VALID replicates.
#
# WHY. Web Table S1 (tab:a2ratio) and every main-text/appendix sentence that states the SIZE of the
# Approach-2 weight-estimation correction were built by dev/a2_ratio_table.R from metrics_ci.csv, i.e.
# from each method's mean SE over ITS OWN valid replicates. M6 and M6-A2 fail on different replicates
# (construction failure 0.348% vs 0.378%; worst cell 812 vs 793 of 1000), so the unpaired % change mixes
# the correction with a change in which datasets were analysed. The paper therefore states the
# size of the correction on the replicates where BOTH SEs are valid and reports failure rates separately.
#
# WHAT. Per (scenario, horizon): mean SE of M6 and of M6-A2 over the replicates where both rows exist and
# both SEs pass the canonical guard of simulation/metrics.R
# (.valid_se: finite, > 0, <= 0.5 for the MRD); % change = 100 * (mean A2 - mean A1) / mean A1.
# Done for the analytic subject-clustered SE (the table body and every pooled figure) and for the
# multiplier SE (the tau=1 resampling-noise statement).
#
# CHECKS FIRST. Section A reproduces the CURRENT (unpaired) Web Table S1 and pooled figures exactly from
# metrics_ci.csv; section B rebuilds the unpaired mean SEs from replicates_all.rds and checks them against
# metrics_ci.csv. Only then is the paired table (section C) trusted.
#
# Sections F and G are side checks: F shows that the -3.6% pooled minimum arises only when the SE <= 0.5
# check is dropped; G re-scores the tab:e2 coverage move ('0.1--0.6 percentage points') on common replicates.
#
# Canonical exact results only. Needs data/exact/results/replicates_all.rds (91 MB; not in git -- see
# the README, "Large files"). Measured runtime 42-61 s on a laptop (about 8 s of it reading that file). Writes nothing in the repo (results RDS to $A2_OUT or tempdir()).
#   Rscript dev/approach2_paired.R
suppressWarnings(suppressMessages({library(stats)}))
ROOT <- Sys.getenv("REPO_ROOT", ".")
RES  <- file.path(ROOT, "data/exact/results")
SE_MAX_MRD <- 0.5
valid_se <- function(se) is.finite(se) & se > 0 & se <= SE_MAX_MRD   # = metrics.R .valid_se("mrd", .)
n_of <- function(s) as.integer(sub("^n(\\d+)_.*", "\\1", s))
fmt2 <- function(x) sprintf("%.2f", x)
stopifnot_eq <- function(got, want, what, tol = 0) {
  ok <- isTRUE(all.equal(got, want, tolerance = tol, check.attributes = FALSE)) ||
        (tol == 0 && identical(unname(got), unname(want)))
  cat(sprintf("  [%s] %s\n", if (ok) "OK" else "MISMATCH", what))
  if (!ok) { print(got); print(want); stop("reproduction check failed: ", what) }
}

summarise_cells <- function(d, label) {
  # d: scenario, n, tau, pct (one row per scenario x horizon)
  p <- d$pct
  cat(sprintf("\n--- %s: %d cells ---\n", label, nrow(d)))
  cat(sprintf("  signed median %+.3f%% | median |change| %.3f%% | range %+.2f%% to %+.2f%%\n",
              median(p), median(abs(p)), min(p), max(p)))
  cat(sprintf("  |change| > 1%%: %d of %d (%.1f%%) | > 2%%: %.1f%% | negative: %d (%.1f%%) | exactly 0: %d\n",
              sum(abs(p) > 1), length(p), 100 * mean(abs(p) > 1), 100 * mean(abs(p) > 2),
              sum(p < 0), 100 * mean(p < 0), sum(p == 0)))
  w <- d[which.min(p), ]; W <- d[which.max(p), ]
  cat(sprintf("  min at %s tau=%g (%+.3f%%); max at %s tau=%g (%+.3f%%)\n",
              w$scenario, w$tau, w$pct, W$scenario, W$tau, W$pct))
  invisible(NULL)
}

grid_tab <- function(d, f = median) {
  t <- tapply(d$pct, list(n = d$n, tau = d$tau), f)
  t[order(as.integer(rownames(t))), order(as.numeric(colnames(t))), drop = FALSE]
}

by_tau <- function(d) {
  do.call(rbind, lapply(split(d, d$tau), function(g) data.frame(
    tau = g$tau[1], median_pct = round(median(g$pct), 2), median_abs = round(median(abs(g$pct)), 2),
    max_abs = round(max(abs(g$pct)), 2), frac_gt1 = round(mean(abs(g$pct) > 1), 3),
    frac_neg = round(mean(g$pct < 0), 3))))
}

# =====================================================================================================
cat("######## A. CURRENT (UNPAIRED) WEB TABLE S1 FROM metrics_ci.csv ########\n")
ci <- read.csv(file.path(RES, "metrics_ci.csv"), stringsAsFactors = FALSE)
ci <- ci[ci$estimand == "mrd", ]
unpaired <- function(m1, m2) {
  a <- ci[ci$method == m1, c("scenario", "tau", "mean_se", "n_sim")]
  b <- ci[ci$method == m2, c("scenario", "tau", "mean_se", "n_sim")]
  d <- merge(a, b, by = c("scenario", "tau"), suffixes = c("1", "2"))
  d <- d[is.finite(d$mean_se1) & is.finite(d$mean_se2) & d$mean_se1 > 0, ]
  d$n <- n_of(d$scenario)
  d$pct <- 100 * (d$mean_se2 - d$mean_se1) / d$mean_se1
  d[order(d$n, d$scenario, d$tau), ]
}
U_an <- unpaired("m6_analytic_cluster", "m6a2_analytic_cluster")
U_mu <- unpaired("m6_multiplier", "m6a2_multiplier")
tabU <- grid_tab(U_an)
print(round(tabU, 2))
# The body of Web Table S1 as first published, from the UNPAIRED means (dev/a2_ratio_table.R).
published <- rbind(c(0.00, -0.52, -1.71, -1.99, -0.71),
                   c(0.00, -0.28, -1.08, -1.66, -1.49),
                   c(0.00, -0.09, -0.51, -0.92, -0.99))
stopifnot_eq(round(tabU, 2) + 0, published, "unpaired table body == published Web Table S1", tol = 1e-12)
stopifnot_eq(nrow(U_an), 405L, "405 analytic cells")
stopifnot_eq(round(median(U_an$pct), 2), -0.35, "pooled signed median -0.35%")
stopifnot_eq(round(c(min(U_an$pct), max(U_an$pct)), 1), c(-3.3, 4.7), "range -3.3% to +4.7%")
stopifnot_eq(round(100 * mean(abs(U_an$pct) > 1)), 36, "36% of cells exceed 1%")
stopifnot_eq(round(100 * mean(U_an$pct < 0)), 77, "77% negative")
stopifnot_eq(round(median(abs(U_an$pct)), 2), 0.42, "median |change| 0.42%")
u1 <- U_an[U_an$tau == 1, ]
stopifnot_eq(sum(u1$pct == 0), 79L, "tau=1 analytic exactly zero in 79 of 81 (unpaired)")
cat(sprintf("  tau=1 analytic non-zero cells (unpaired): %s\n",
            paste(sprintf("%s %+.4f%% (n_sim %d vs %d)", u1$scenario[u1$pct != 0], u1$pct[u1$pct != 0],
                          u1$n_sim1[u1$pct != 0], u1$n_sim2[u1$pct != 0]), collapse = "; ")))
m1 <- U_mu[U_mu$tau == 1, ]
cat(sprintf("  tau=1 multiplier (unpaired): exactly zero in %d of %d; max |change| %.3f%%\n",
            sum(m1$pct == 0), nrow(m1), max(abs(m1$pct))))
summarise_cells(U_an, "UNPAIRED analytic (as published)")
summarise_cells(U_mu, "UNPAIRED multiplier")
cat("\nunpaired median |change| by n x tau (analytic):\n"); print(round(grid_tab(U_an, function(x) median(abs(x))), 2))

# =====================================================================================================
cat("\n######## B. REBUILD FROM replicates_all.rds ########\n")
t0 <- Sys.time()
r <- readRDS(file.path(RES, "replicates_all.rds"))
METH <- c("m6_analytic_cluster", "m6a2_analytic_cluster", "m6_multiplier", "m6a2_multiplier")
r <- r[r$estimand == "mrd" & r$method %in% METH,
       c("scenario", "n", "tau", "rep", "method", "se", "lower", "upper")]
cat(sprintf("  read + filter: %.1f s; %d rows\n", as.numeric(Sys.time() - t0, units = "secs"), nrow(r)))
r$key <- paste(r$scenario, r$tau, r$rep, sep = "|")
for (m in METH) stopifnot(!anyDuplicated(r$key[r$method == m]))
cat("  [OK] one row per (scenario, tau, rep) per method\n")
r$ok <- valid_se(r$se)

# B1: unpaired mean SE from the replicates must reproduce metrics_ci.csv mean_se
own <- aggregate(se ~ scenario + tau + method, data = r[r$ok, ], FUN = mean)
chk <- merge(own, ci[ci$method %in% METH, c("scenario", "tau", "method", "mean_se")],
             by = c("scenario", "tau", "method"))
stopifnot_eq(nrow(chk), 4L * 405L, "1620 (scenario, tau, method) cells matched")
stopifnot_eq(max(abs(chk$se - chk$mean_se)), 0, "replicate means == metrics_ci mean_se", tol = 1e-12)
cat(sprintf("  invalid-but-present SE rows (guard): %s\n",
            paste(sprintf("%s=%d", METH, sapply(METH, function(m) sum(!r$ok[r$method == m]))),
                  collapse = " ")))

# =====================================================================================================
cat("\n######## C. PAIRED: COMMON VALID REPLICATES ########\n")
pair_reps <- function(m1, m2) {
  a <- r[r$method == m1 & r$ok, c("key", "scenario", "n", "tau", "rep", "se")]
  b <- r[r$method == m2 & r$ok, c("key", "se")]
  merge(a, b, by = "key", suffixes = c("1", "2"))
}
paired_cells <- function(p) {
  g <- split(p, paste(p$scenario, p$tau))
  d <- do.call(rbind, lapply(g, function(x) data.frame(
    scenario = x$scenario[1], n = x$n[1], tau = x$tau[1], n_common = nrow(x),
    se1 = mean(x$se1), se2 = mean(x$se2), identical_reps = sum(x$se1 == x$se2),
    stringsAsFactors = FALSE)))
  d$pct <- 100 * (d$se2 - d$se1) / d$se1
  rownames(d) <- NULL
  d[order(d$n, d$scenario, d$tau), ]
}
Pr_an <- pair_reps("m6_analytic_cluster", "m6a2_analytic_cluster")
Pr_mu <- pair_reps("m6_multiplier", "m6a2_multiplier")
P_an <- paired_cells(Pr_an)
P_mu <- paired_cells(Pr_mu)
stopifnot_eq(nrow(P_an), 405L, "405 paired analytic cells")
stopifnot_eq(nrow(P_mu), 405L, "405 paired multiplier cells")

# how much pairing drops, per method
dropinfo <- function(P, m1, m2) {
  n1 <- tapply(r$ok[r$method == m1], paste(r$scenario, r$tau)[r$method == m1], sum)
  n2 <- tapply(r$ok[r$method == m2], paste(r$scenario, r$tau)[r$method == m2], sum)
  k <- paste(P$scenario, P$tau)
  cat(sprintf("  %s vs %s: common replicate-horizons %d; A1-only valid %d, A2-only valid %d; cells where the sets differ: %d of %d\n",
              m1, m2, sum(P$n_common), sum(n1[k]) - sum(P$n_common), sum(n2[k]) - sum(P$n_common),
              sum(n1[k] != P$n_common | n2[k] != P$n_common), nrow(P)))
  worst <- which.max(pmax(n1[k], n2[k]) - P$n_common)
  cat(sprintf("    largest set difference: %s tau=%g: A1 %d, A2 %d, common %d\n",
              P$scenario[worst], P$tau[worst], n1[k][worst], n2[k][worst], P$n_common[worst]))
}
dropinfo(P_an, "m6_analytic_cluster", "m6a2_analytic_cluster")
dropinfo(P_mu, "m6_multiplier", "m6a2_multiplier")

cat("\nPAIRED analytic: median % change by n x tau (THE NEW TABLE BODY)\n")
tabP <- grid_tab(P_an); print(round(tabP, 2))
cat("\nchange vs published body (paired - unpaired, rounded):\n"); print(round(tabP, 2) - round(tabU, 2))
cat("\nPAIRED analytic: median |change| by n x tau\n"); print(round(grid_tab(P_an, function(x) median(abs(x))), 2))
cat("\nPAIRED analytic: by horizon (pooled over n)\n"); print(by_tau(P_an))
cat("\nPAIRED analytic: by n (pooled over tau)\n")
print(do.call(rbind, lapply(split(P_an, P_an$n), function(g) data.frame(
  n = g$n[1], median_pct = round(median(g$pct), 3), max_abs = round(max(abs(g$pct)), 2),
  frac_gt1 = round(mean(abs(g$pct) > 1), 3)))))
summarise_cells(P_an, "PAIRED analytic (common valid replicates)")
summarise_cells(P_mu, "PAIRED multiplier (common valid replicates)")
cat("\nPAIRED multiplier: median % change by n x tau\n"); print(round(grid_tab(P_mu), 2))

cat("\n--- tau = 1 under pairing ---\n")
a1 <- P_an[P_an$tau == 1, ]; x1 <- P_mu[P_mu$tau == 1, ]
cat(sprintf("  analytic : %% change exactly 0 in %d of %d scenarios; max |change| %.2e%%\n",
            sum(a1$pct == 0), nrow(a1), max(abs(a1$pct))))
cat(sprintf("             per-replicate SEs bitwise identical in %d of %d replicate pairs (%d scenarios all-identical)\n",
            sum(a1$identical_reps), sum(a1$n_common), sum(a1$identical_reps == a1$n_common)))
cat(sprintf("  multiplier: %% change exactly 0 in %d of %d scenarios; max |change| %.3f%%; median |change| %.3f%%\n",
            sum(x1$pct == 0), nrow(x1), max(abs(x1$pct)), median(abs(x1$pct))))

# =====================================================================================================
cat("\n######## D. SENTENCE CHECKS (appendix B growth sentence) ########\n")
mag <- grid_tab(P_an, function(x) median(abs(x)))
sgn <- tabP
for (nn in rownames(mag)) {
  cat(sprintf("  n=%s  median |change| tau1..5: %s | signed: %s\n", nn,
              paste(fmt2(mag[nn, ]), collapse = " "), paste(fmt2(sgn[nn, ]), collapse = " ")))
}
grows_to4 <- all(apply(abs(sgn[, 1:4]), 1, function(v) all(diff(v) > 0)))
cat(sprintf("  |signed median| strictly increases tau=1..4 at every n: %s\n", grows_to4))
cat(sprintf("  |signed median| at tau=5 vs tau=4: %s\n",
            paste(sprintf("n=%s %s", rownames(sgn), ifelse(abs(sgn[, 5]) > abs(sgn[, 4]), "grows", "falls back")),
                  collapse = "; ")))
cat(sprintf("  median |change| at tau=5 vs tau=4: %s\n",
            paste(sprintf("n=%s %s->%s", rownames(mag), fmt2(mag[, 4]), fmt2(mag[, 5])), collapse = "; ")))
cat("  dispersion across the 27 scenarios (IQR of % change) by n x tau:\n")
print(round(grid_tab(P_an, IQR), 2))
cat("  dispersion across the 27 scenarios (SD of % change) by n x tau:\n")
print(round(grid_tab(P_an, sd), 2))
cat(sprintf("  cells with |change|>4%%: %s\n",
            paste(sprintf("%s tau=%g %+.2f", P_an$scenario, P_an$tau, P_an$pct)[abs(P_an$pct) > 4], collapse = "; ")))
cat(sprintf("  tau of the 10 largest |change|: %s\n",
            paste(P_an$tau[order(-abs(P_an$pct))][1:10], collapse = " ")))

# =====================================================================================================
cat("\n######## E. LaTeX: regenerated Web Table S1 body ########\n")
for (nn in rownames(tabP)) {
  cells <- sapply(tabP[nn, ], function(v) if (round(v, 2) == 0) "$0.00$" else sprintf("$%.2f$", v))
  cat(sprintf("$%s$%s & %s\\\\\n", nn, strrep(" ", 6 - nchar(nn)), paste(cells, collapse = " & ")))
}
cat(sprintf("\nPOOLED (paired): median %.2f%%, range %.1f%% to %+.1f%%, %.0f%% exceed 1%%, %.0f%% negative, median |change| %.2f%%\n",
            median(P_an$pct), min(P_an$pct), max(P_an$pct), 100 * mean(abs(P_an$pct) > 1),
            100 * mean(P_an$pct < 0), median(abs(P_an$pct))))

# =====================================================================================================
cat("\n######## F. WHY NOT -3.6%: the same pairing WITHOUT the SE <= 0.5 check ########\n")
# The paper's mean estimated SE excludes SEs above 0.5 (manuscript 'Validity checks on scale statistics';
# Web Appendix D). Dropping that check (finite and > 0 only) changes only the pooled minimum.
r$fin <- is.finite(r$se) & r$se > 0
ug <- merge(r[r$method == "m6_analytic_cluster" & r$fin, c("key", "scenario", "n", "tau", "se")],
            r[r$method == "m6a2_analytic_cluster" & r$fin, c("key", "se")], by = "key")
ug <- aggregate(cbind(se.x, se.y) ~ scenario + n + tau, data = ug, FUN = mean)
ug$pct <- 100 * (ug$se.y - ug$se.x) / ug$se.x
summarise_cells(ug, "PAIRED analytic, finite & > 0 only (no SE <= 0.5 check)")
cat("  table body unchanged by the check: ",
    isTRUE(all.equal(round(grid_tab(ug), 2), round(tabP, 2))), "\n")

# =====================================================================================================
cat("\n######## G. COVERAGE MOVE (tab:e2 '0.1--0.6 percentage points') on common replicates ########\n")
em <- read.csv(file.path(RES, "metrics_empirical.csv"), stringsAsFactors = FALSE)
em <- em[em$estimand == "mrd", ]
r$truth <- em$truth[match(paste(r$scenario, r$tau), paste(em$scenario, em$tau))]
r$cov <- r$lower <= r$truth & r$upper >= r$truth
for (pr in list(METH[1:2], METH[3:4])) {
  a <- r[r$method == pr[1] & is.finite(r$cov), c("key", "scenario", "n", "tau", "cov")]
  b <- r[r$method == pr[2] & is.finite(r$cov), c("key", "cov")]
  d <- merge(a, b, by = "key")
  cell <- aggregate(cbind(cov.x, cov.y) ~ scenario + n + tau, data = d, FUN = mean)
  byn <- aggregate(cbind(cov.x, cov.y) ~ n, data = cell, FUN = mean)
  cat(sprintf("  %-22s -> %-24s A1-A2 (pp) at n=300/1000/5000: %s\n", pr[1], pr[2],
              paste(sprintf("%.2f", 100 * (byn$cov.x - byn$cov.y)), collapse = " / ")))
}

saveRDS(list(unpaired_analytic = U_an, unpaired_multiplier = U_mu,
             paired_analytic = P_an, paired_multiplier = P_mu, table_body = tabP),
        file.path(Sys.getenv("A2_OUT", tempdir()), "approach2_paired_results.rds"))
cat("\nAPPROACH-2 PAIRED DONE\n")
