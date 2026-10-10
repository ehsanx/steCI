# ############################################################################
# ## ⚠️  THIS SCRIPT'S HEADLINE CONCLUSION WAS **REFUTED**. READ THIS FIRST. ##
# ############################################################################
#
# WHAT THIS SCRIPT GOT WRONG:
#  * "ESS collapses to 14 of 930" is a **23-FOLD ERROR** from a 3-SEED LOW-TAIL DRAW. With 24 seeds
#    the same cell gives ess_frac = 0.326 (MC SE 0.055); per-seed range 0.0026-0.7895 (307-fold);
#    the reported 0.015 is the 4th smallest of 24 (P(<=0.05) = 0.022). It does not resolve at ANY
#    seed count run so far (8 fresh seeds -> 0.119 vs 24 seeds -> 0.326). **Report no ESS value.**
#  * The 1.000 baseline is an **algebraic identity**, not a measurement: bin (0,1] is exclusively
#    rownum==1, where ste_prepare.R:81-82 sets ipw.stab.cum = 1 by construction.
#  * The "horizon bins" are NOT horizon bins. Every risk set lies in exactly one rownum stratum
#    (rownum == ceiling(t)), so the bins index HOW MANY RATIO FACTORS have been multiplied (0..4).
#  * ESS is NOT the mechanism, and anti-predicts the shortfall: holding gamma.L=3 fixed, a 7.2-fold
#    ESS swing moves se_ratio_fold 1.09-fold IN THE WRONG DIRECTION (within-level r = -0.793).
#    cor(gamma.L, se_ratio_fold) = -0.948 vs cor(ess_frac, se_ratio_fold) = 0.534 -- ESS is a
#    redundant, WORSE relabelling of gamma.L.
#  * The saturation framing is a non sequitur: Kish ESS is a function of w ALONE -- Dinv is not in it.
#    Dinv is the MITIGATION (a ceiling on max_share), not the cause. And ci_lef.R:98 makes the SE the
#    root-sum-square of per-subject IFs, so within-sample concentration is ALREADY in the SE: a
#    14-subject average implies a LARGE SE and CONSERVATIVE coverage. Ours is too SMALL. The missing
#    variance is BETWEEN-sample; ESS is WITHIN-sample.
#  * E[W^2] < infinity (proved separately) implies ess_frac -> (EW)^2/E[W^2] > 0, a positive constant,
#    so a collapse to 0.015-and-falling contradicts our own algebra.
#  * max_share n-growth (3.24x) is an ARTEFACT: 40 seeds give 0.743, CI [0.454, 1.156] (shrinkage).
#
# WHAT SURVIVES AND IS STILL WORTH RUNNING THIS FOR:
#  * C4-3: **max cond(X'WX) = 78.7** across all event times/scenarios -- benign (ill-conditioning is
#    ~1e8). The earlier "D(t)^-1 blows up as the risk set shrinks" hypothesis is DEAD ON MEASUREMENT.
#  * C4-4 read ONLY as a seed-noisy description of the WEIGHT REGIME (control arm only, no point
#    value, explicit MC error) -- never as an explanation of the shortfall.
#
# BEFORE TRUSTING ANY FUTURE NUMBER FROM A SCRIPT LIKE THIS: >= 24 seeds, and clear the four-part
# reliability checklist.
# ############################################################################
#
# The Dinv self-normalisation question. Cheap, direct, no re-run.
#
# THE REFEREE QUESTION we currently cannot answer: "your estimator self-normalises through
# Dinv = (X'WX)^-1 (aalen_if.R:47), so subject i's contribution behaves like
# w_i / (sum_j w_j + w_i), which SATURATES at 1 and is therefore BOUNDED. How can heavy weights
# matter at all?"  One objection is that this is a contradiction: our mechanism story and the
# self-normalisation explanation cannot both be true.
#
# THIS IS DECIDABLE BY MEASUREMENT, not argument. Saturation only bites if some subject is ACTUALLY
# near the saturation regime, i.e. if its weight share of the risk set is O(1). So measure:
#
#   (1) SHARE  = max_i w_i / sum_{j in risk set} w_j, by event time and by n.
#       -> If max share is small (<<1) NO subject is anywhere near saturation, the objection
#          dissolves, and the effect must come from AGGREGATE weight dispersion, not from any one
#          subject dominating. If max share is O(1), saturation is real and the story needs restating.
#       -> Its scaling in n is the discriminator: a share ~ 1/n means ordinary averaging; a share
#          that does NOT shrink like 1/n means a few subjects carry the fit.
#
#   (2) COND(X'WX) by event time / risk-set size -- the ORIGINAL "D(t)^-1 is unstable as the
#       risk set shrinks" hypothesis, which was argued about earlier but NEVER DIRECTLY MEASURED.
#
#   (3) ESS of the risk set = (sum w)^2 / sum w^2, by event time -- how many subjects the weighted
#       fit is effectively using at late times.
#
#   Rscript dev/phaseC_dinv.R
suppressWarnings(suppressMessages({library(dplyr)}))
ROOT <- "."
for (f in list.files(file.path(ROOT, "steCI/R"), pattern = "[.]R$", full.names = TRUE)) source(f)

# Recompute the per-event-time quantities that .aalen_if() forms internally at aalen_if.R:45-47.
probe <- function(n, gamma.L, gamma.0, alpha.0, seed) {
  d  <- dat_sim(n = n, n.visit = 5, gamma.0 = gamma.0, gamma.L = gamma.L, alpha.0 = alpha.0,
                alpha.A = -0.04, alpha.L = 0.015, alpha.U = 0.015, seed = seed)
  ds <- ste_prepare(d, n.visit = 5, truncate = 0.95)
  X <- cbind(1, ds$A.baseline, ds$L.baseline)
  w <- ds$ipw.stab.cum
  start <- ds$time.new; stop <- ds$time.stop.new; ev <- ds$event == 1
  etimes <- sort(unique(stop[ev]))
  out <- lapply(etimes, function(tj) {
    inR <- start < tj & stop >= tj
    if (sum(inR) < 5) return(NULL)
    Xr <- X[inR, , drop = FALSE]; wr <- w[inR]
    sw <- sum(wr)
    XWX <- crossprod(Xr, wr * Xr)
    cn <- tryCatch(kappa(XWX, exact = TRUE), error = function(e) NA_real_)
    data.frame(t = tj, n_risk = sum(inR),
               max_share = max(wr) / sw,                    # (1) the saturation diagnostic
               share_x_nrisk = max(wr) / sw * sum(inR),     #     = 1 under equal weights
               ess_risk = sw^2 / sum(wr^2),                 # (3)
               ess_frac = (sw^2 / sum(wr^2)) / sum(inR),
               cond_XWX = cn,                               # (2)
               max_w = max(wr))
  })
  out <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
  out$n <- n; out$gamma.L <- gamma.L; out
}

cat("######## does the fit self-normalise into SATURATION? ########\n\n")
cat("Reading aalen_if.R:45-53 -- at each event time tj the row contribution is\n")
cat("   IF_r += Dinv %*% X_r * w_r * dM_r ,  Dinv = (X_R' W_R X_R)^-1 on the risk set R.\n")
cat("So subject r's weight enters the numerator AND the denominator sum => the share\n")
cat("   w_r / sum_{j in R} w_j   saturates at 1 IF some w_r dominates. Does one?\n\n")

grid <- expand.grid(n = c(300L, 1000L, 5000L), gamma.L = c(0.5, 3.0), seed = 1:3)
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
  with(grid[i, ], probe(n, gamma.L, -2.0, 0.25, seed))))

cat("C4-1. THE SATURATION DIAGNOSTIC: max_i w_i / sum_j w_j on the risk set\n")
cat("      (saturation is only a live concern if this is O(1). Under EQUAL weights it is 1/n_risk.)\n")
s1 <- res %>% group_by(gamma.L, n) %>%
  summarise(median_max_share = round(median(max_share), 5),
            worst_max_share  = round(max(max_share), 4),
            median_n_risk    = round(median(n_risk)),
            `worst_share_x_n_risk` = round(max(share_x_nrisk), 2), .groups = "drop")
print(as.data.frame(s1))
cat("\n   Read: share_x_n_risk = (max share) x (risk-set size) = 1.0 under equal weights.\n")
cat("   A value of k means the heaviest subject carries k subjects' worth of the fit.\n")
cat("   For SATURATION to be the story, max_share must approach 1 -- i.e. share_x_n_risk ~ n_risk.\n")

cat("\nC4-2. does max_share shrink like 1/n?  (ordinary averaging => yes)\n")
s2 <- res %>% group_by(gamma.L, n) %>%
  summarise(mean_max_share = mean(max_share), .groups = "drop") %>%
  group_by(gamma.L) %>% arrange(n, .by_group = TRUE) %>%
  mutate(ratio_vs_prev = round(mean_max_share / lag(mean_max_share), 3),
         mean_max_share = signif(mean_max_share, 4))
print(as.data.frame(s2))
cat("   n goes 300 -> 1000 (x3.33) -> 5000 (x5). Exact 1/n scaling => ratios 0.30 and 0.20.\n")

cat("\nC4-3. THE D(t)^-1 HYPOTHESIS, MEASURED AT LAST: cond(X'WX) vs risk-set size\n")
cat("      (earlier we conjectured D(t)^-1 blows up as the risk set shrinks; argued but never measured.)\n")
res$t_bin <- cut(res$t, breaks = c(0, 1, 2, 3, 4, 5), labels = c("(0,1]", "(1,2]", "(2,3]", "(3,4]", "(4,5]"))
s3 <- res %>% filter(!is.na(t_bin)) %>% group_by(gamma.L, t_bin) %>%
  summarise(median_n_risk = round(median(n_risk)),
            median_cond = round(median(cond_XWX, na.rm = TRUE), 1),
            p95_cond = round(quantile(cond_XWX, 0.95, na.rm = TRUE), 1),
            max_cond = round(max(cond_XWX, na.rm = TRUE), 1), .groups = "drop")
print(as.data.frame(s3))
cat("   Read: a condition number in the tens is benign; ill-conditioning means 1e8+.\n")

cat("\nC4-4. effective sample size of the RISK SET by horizon (how much data the fit really has)\n")
s4 <- res %>% filter(!is.na(t_bin)) %>% group_by(gamma.L, t_bin) %>%
  summarise(median_n_risk = round(median(n_risk)),
            median_ess = round(median(ess_risk)),
            median_ess_frac = round(median(ess_frac), 3), .groups = "drop")
print(as.data.frame(s4))

saveRDS(res, file.path(ROOT, "dev/phaseC_dinv_results.rds"))
cat("\nINTERPRETATION KEY:\n")
cat(" * max_share << 1 and shrinking like 1/n  -> NO subject is near saturation. The referee's\n")
cat("   objection dissolves: self-normalisation WOULD bound a dominant subject, but none exists.\n")
cat("   The effect then comes from AGGREGATE weight dispersion, not from one subject dominating,\n")
cat("   and the mechanism story and the self-normalisation are CONSISTENT (not contradictory).\n")
cat(" * max_share ~ O(1)  -> saturation is real and the mechanism story must be restated.\n")
cat(" * cond(X'WX) benign at late t -> the earlier D(t)^-1 hypothesis is dead on measurement too.\n")
cat("\nSELF-NORMALISATION PROBE DONE (saved dev/phaseC_dinv_results.rds)\n")
