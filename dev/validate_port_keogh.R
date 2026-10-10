# PORT VALIDATION: does steCI::ste_prepare() reproduce Keogh's analysis_seq.R weights EXACTLY?
#
# WHY THIS EXISTS. The earlier result (the estimator is not root-n over n in [300,5000]; the SE
# shortfall grows with horizon) rests entirely on ipw.stab.cum being CORRECT. ste_prepare() is a
# clean-room port of Keogh et al. (2023) `sequential_trials/analysis_seq.R`. If the port makes the
# weights more extreme than the reference, then "the estimator is not root-n" is really "our port is
# broken" and that result evaporates. This script tests the port against the reference NUMERICALLY.
#
# METHOD. Generate dat.long once with steCI::dat_sim(), then build the weights TWICE:
#   (a) via steCI::ste_prepare()
#   (b) via a VERBATIM transcription of analysis_seq.R lines 69-182
# and compare ipw.stab.cum element-by-element on the (id, trial, rownum) key.
#
# DEVIATIONS FROM VERBATIM (both semantically inert; documented so the test is honest):
#   1. rbind() -> dplyr::bind_rows(). Keogh's trial.0 chain never creates Alag1.baseline (only
#      trials 1-4 need it, for the eligibility filter), so rbind() would fail on mismatched columns.
#      bind_rows fills NA. Alag1.baseline is never referenced after the per-trial filter, so this
#      cannot affect any weight.
#   2. as.data.frame() before the weight block, because Keogh's `dat.seq[cond,]$col <- ...`
#      subassignment idiom does not work on a grouped tibble. ste_prepare() does the same at :54.
# Everything else -- the glm formulas, the `rownum != 5` predrows rule, the A.baseline==1 override,
# the lag/cumprod, the 0.95 truncation -- is transcribed character-for-character.
#
#   Rscript dev/validate_port_keogh.R
suppressWarnings(suppressMessages({library(dplyr); library(tidyr)}))
ROOT <- "."
for (f in list.files(file.path(ROOT, "steCI/R"), pattern = "[.]R$", full.names = TRUE)) source(f)

# ---------------------------------------------------------------- Keogh, verbatim
keogh_weights <- function(dat.long, n.visit = 5) {
  mk <- function(g) {
    if (g == 0) {
      trial <- dat.long
      trial <- trial %>%
        mutate(trial = 0) %>%
        group_by(id) %>%
        mutate(ones = 1) %>% mutate(rownum = cumsum(ones)) %>%
        mutate(A.baseline = first(A)) %>%
        mutate(L.baseline = first(L)) %>%
        mutate(time.new = time - trial) %>%
        mutate(time.stop.new = time.stop - trial)
    } else {
      trial <- dat.long[dat.long$time >= g, ]
      trial <- trial %>%
        mutate(trial = g) %>%
        group_by(id) %>%
        mutate(ones = 1) %>% mutate(rownum = cumsum(ones)) %>%
        mutate(time.new = time - trial) %>%
        mutate(time.stop.new = time.stop - trial) %>%
        mutate(Alag1.baseline = first(Alag1)) %>%
        filter(Alag1.baseline == 0) %>%   # restrict to those not previously treated at trial start
        mutate(A.baseline = first(A)) %>%
        mutate(L.baseline = first(L))
    }
    trial
  }
  dat.seq.orig <- bind_rows(lapply(0:(n.visit - 1), mk))   # DEVIATION 1 (see header)

  dat.seq.orig <- dat.seq.orig %>% group_by(id, trial) %>% mutate(A.lead1 = lead(A, n = 1))
  dat.seq.orig <- dat.seq.orig %>% group_by(id, trial) %>% mutate(L.lead1 = lead(L, n = 1))
  dat.seq.orig <- dat.seq.orig %>% mutate(Anext.equal.to.baseline = ifelse(A.lead1 == A.baseline, 1, 0))

  dat.seq <- dat.seq.orig[dat.seq.orig$A == dat.seq.orig$A.baseline, ]
  dat.seq <- as.data.frame(dat.seq)                        # DEVIATION 2 (see header)

  wt.mod.denom <- glm(Anext.equal.to.baseline ~ L.lead1, family = "binomial",
                      data = dat.seq[dat.seq$A.baseline == 0, ])
  dat.seq$wt.denom <- 1
  dat.seq[dat.seq$rownum != 5, ]$wt.denom <-
    predict(wt.mod.denom, newdata = dat.seq[dat.seq$rownum != 5, ], type = "response")
  dat.seq$wt.denom <- ifelse(dat.seq$A.baseline == 1, 1, dat.seq$wt.denom)

  wt.mod.num <- glm(Anext.equal.to.baseline ~ as.factor(rownum) * L.baseline, family = "binomial",
                    data = dat.seq[dat.seq$A.baseline == 0, ])
  dat.seq$wt.num <- 1
  dat.seq[dat.seq$rownum != 5, ]$wt.num <-
    predict(wt.mod.num, newdata = dat.seq[dat.seq$rownum != 5, ], type = "response")
  dat.seq$wt.num <- ifelse(dat.seq$A.baseline == 1, 1, dat.seq$wt.num)

  dat.seq$ipw.stab <- dat.seq$wt.num / dat.seq$wt.denom

  dat.seq <- dat.seq %>% group_by(id, trial) %>% mutate(ipw.stab.lag = lag(ipw.stab, n = 1))
  dat.seq <- dat.seq %>% mutate(ipw.stab.lag = ifelse(rownum == 1, 1, ipw.stab.lag))
  dat.seq <- dat.seq %>% group_by(id, trial) %>% mutate(ipw.stab.cum = cumprod(ipw.stab.lag))

  pct95 <- quantile(dat.seq$ipw.stab.cum, 0.95)
  dat.seq$ipw.stab.cum.trunc <- ifelse(dat.seq$ipw.stab.cum > pct95, pct95, dat.seq$ipw.stab.cum)
  as.data.frame(dat.seq)
}

# ---------------------------------------------------------------- compare on one dataset
compare_one <- function(n, alpha.0, gamma.0, gamma.L, seed) {
  d <- dat_sim(n = n, n.visit = 5, gamma.0 = gamma.0, gamma.L = gamma.L,
               alpha.0 = alpha.0, alpha.A = -0.04, alpha.L = 0.015, alpha.U = 0.015, seed = seed)
  ours  <- ste_prepare(d, n.visit = 5, truncate = 0.95)
  theirs <- keogh_weights(d, n.visit = 5)

  k <- function(x) paste(x$id, x$trial, x$rownum, sep = "|")
  ko <- k(ours); kt <- k(theirs)
  if (anyDuplicated(ko) || anyDuplicated(kt)) stop("non-unique key")

  m <- match(ko, kt)
  cmp <- function(a, b) { v <- abs(a - b); max(v[is.finite(v)], na.rm = TRUE) }
  na_mismatch <- function(a, b) sum(is.na(a) != is.na(b))

  data.frame(
    n = n, alpha.0 = alpha.0, gamma.0 = gamma.0, gamma.L = gamma.L, seed = seed,
    rows_ours = nrow(ours), rows_theirs = nrow(theirs),
    rows_unmatched = sum(is.na(m)),
    d_wt_denom  = cmp(ours$wt.denom,          theirs$wt.denom[m]),
    d_wt_num    = cmp(ours$wt.num,            theirs$wt.num[m]),
    d_ipw_stab  = cmp(ours$ipw.stab,          theirs$ipw.stab[m]),
    d_CUM       = cmp(ours$ipw.stab.cum,      theirs$ipw.stab.cum[m]),
    d_CUM_TRUNC = cmp(ours$ipw.stab.cum.trunc, theirs$ipw.stab.cum.trunc[m]),
    na_cum      = na_mismatch(ours$ipw.stab.cum, theirs$ipw.stab.cum[m]),
    max_w_ours  = max(ours$ipw.stab.cum, na.rm = TRUE),
    max_w_theirs = max(theirs$ipw.stab.cum, na.rm = TRUE),
    stringsAsFactors = FALSE)
}

cat("######## PORT VALIDATION: ste_prepare() vs Keogh analysis_seq.R ########\n\n")
cat("Comparing ipw.stab.cum element-by-element on the (id, trial, rownum) key.\n")
cat("A faithful port => every d_* column is 0 (or ~1e-16) and rows_unmatched = 0.\n\n")

grid <- expand.grid(n = c(300L, 1000L), alpha.0 = c(0.15, 0.25, 0.40),
                    gamma.0 = c(0.0, -2.0), gamma.L = c(0.5, 3.0),
                    seed = 1:3, stringsAsFactors = FALSE)
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i)
  with(grid[i, ], compare_one(n, alpha.0, gamma.0, gamma.L, seed))))

cat("--- worst discrepancy across all", nrow(res), "datasets ---\n")
for (cc in c("d_wt_denom", "d_wt_num", "d_ipw_stab", "d_CUM", "d_CUM_TRUNC")) {
  cat(sprintf("  %-12s max = %s\n", cc, format(max(res[[cc]], na.rm = TRUE), digits = 4)))
}
cat("  rows_unmatched total =", sum(res$rows_unmatched), "\n")
cat("  NA-pattern mismatches on ipw.stab.cum =", sum(res$na_cum), "\n")
cat("  row-count mismatches =", sum(res$rows_ours != res$rows_theirs), "\n\n")

cat("--- the hard corner (gamma.L=3, gamma.0=-2), where the finding lives ---\n")
print(as.data.frame(res %>% filter(gamma.L == 3, gamma.0 == -2) %>%
  select(n, alpha.0, seed, d_CUM, d_CUM_TRUNC, max_w_ours, max_w_theirs) %>%
  mutate(across(where(is.numeric), ~signif(., 6)))))

cat("\n--- per-cell detail (all rows) ---\n")
print(as.data.frame(res %>% select(n, alpha.0, gamma.0, gamma.L, seed, rows_ours, rows_theirs,
                                   d_CUM, max_w_ours) %>%
  mutate(across(c(d_CUM, max_w_ours), ~signif(., 6)))))

VERDICT <- if (max(res$d_CUM, na.rm = TRUE) < 1e-10 && sum(res$rows_unmatched) == 0 &&
               sum(res$rows_ours != res$rows_theirs) == 0) "PORT IS FAITHFUL" else "PORT DIVERGES"
cat("\n=========================================================\n")
cat("VERDICT:", VERDICT, "\n")
cat("=========================================================\n")

saveRDS(res, file.path(ROOT, "dev/validate_port_keogh_results.rds"))
cat("\nPORT_VALIDATION DONE (saved dev/validate_port_keogh_results.rds)\n")
