# THE FINITE-MOMENT NEGATIVE CONTROL -- the artifact that was missing.
#
# WHY THIS FILE EXISTS. The tail-index mechanism (candidate #3 of six) was killed by a negative
# control, but the computation was never stored -- it survived only as a one-line assertion in
# a comment in truncation_arc.R. This is BLOCKING: if the control is unrecorded, the
# conclusion rests on nothing. This reproduces it.
#
# THE CONTROL. Take the inverse-propensity primitive of our OWN data-generating mechanism,
#     w_nc = 1 + exp(gamma.0 + gamma.L * L),      L ~ Gaussian (dat_sim.R:43,46)
# and note that EVERY MOMENT OF w_nc IS FINITE BY ALGEBRA: for L ~ N(mu, s^2),
#     E[w_nc^m] = E[(1 + e^{g0 + gL L})^m] = sum_{j=0..m} C(m,j) e^{j*g0} E[e^{j*gL*L}]
#               = sum_{j=0..m} C(m,j) e^{j*g0} exp(j*gL*mu + j^2*gL^2*s^2/2)  <  infinity.
# It is a finite sum of lognormal moments. There is NO tail index: the tail is lognormal-type, so
# alpha = infinity and no power law exists.
#
# THE TEST. Run the standard heavy-tail diagnostics on BOTH the real analysis weight
# (ipw.stab.cum) and w_nc, on the SAME simulated data. If the diagnostics return similar
# "evidence of a heavy tail" for a variable that provably has none, then those diagnostics carry NO
# pathology content and cannot support a tail-index claim about the real weights.
#
# Diagnostics: Hill estimator at several k; mean-excess slope; ESS/n scaling in n; max-weight scaling
# in n. (These are exactly the four the tail-index mechanism relied on.)
#
# RESULT (24 seeds, gamma.L=3, gamma.0=-2, alpha.0=0.25): ALL FOUR DIAGNOSTICS FIRE ON THE CONTROL.
#   Hill (1%/2%/5%): control 0.79/0.73/0.66 vs real 1.48/1.61/2.15   <- control looks WORSE
#   ESS/n slope in n: control -0.521 vs real -0.358                  <- control collapses FASTER
#   max-weight slope: control +0.714 vs real +0.894 ; mean excess rising for both
# A Hill estimate of 0.66-0.80 is deep in the conventional "infinite variance" region (alpha <= 2) for
# a variable whose every moment is finite BY ALGEBRA. Hence a Hill estimate below 2 is NOT evidence of
# infinite variance here, and the tail-index mechanism is dead on evidence rather than assertion.
#
# HONEST CAVEAT: this control is NOT scale-matched. Its log-dispersion is gamma.L*sd(L) ~ 3, against
# sd(log w) ~ 0.31 for the real analysis weight -- roughly a tenfold difference -- and its maximum is
# ~10^5 against ~10^2. The logic is unaffected (exhibiting ANY finite-moment variable with Hill ~ 0.7
# suffices to show that Hill < 2 does not establish infinite variance), but a scale-matched control
# would be a stronger test and is the obvious response if a referee presses.
#
#   Rscript dev/negative_control.R          (NC_SEEDS overrides; default 24, the minimum seed count)
suppressWarnings(suppressMessages({library(dplyr)}))
ROOT <- "."
for (f in list.files(file.path(ROOT, "steCI/R"), pattern = "[.]R$", full.names = TRUE)) source(f)
NSEED <- as.integer(Sys.getenv("NC_SEEDS", "24"))
NS <- c(1000L, 5000L, 20000L)
GL <- 3.0; G0 <- -2.0; A0 <- 0.25

# Hill estimator of the tail index from the top k order statistics.
# alpha_hat = 1 / mean(log x_(i) / x_(k+1)), i = 1..k  (x sorted decreasing)
hill <- function(x, k) {
  x <- sort(x[is.finite(x) & x > 0], decreasing = TRUE)
  if (length(x) <= k + 1) return(NA_real_)
  1 / mean(log(x[1:k] / x[k + 1]))
}
# Mean excess e(u) = E[X - u | X > u]; linear-increasing in u => heavy (Pareto-like) tail.
mean_excess <- function(x, u) { y <- x[x > u]; if (length(y) < 20) NA_real_ else mean(y - u) }
ess_frac <- function(w) sum(w)^2 / (length(w) * sum(w^2))

cat("######## THE FINITE-MOMENT NEGATIVE CONTROL ########\n\n")
cat("Real weight  : ipw.stab.cum from ste_prepare() (cumulative product over visits)\n")
cat("Control      : w_nc = 1 + exp(gamma.0 + gamma.L*L), L Gaussian  => ALL MOMENTS FINITE by algebra\n")
cat("Cells        : gamma.L =", GL, ", gamma.0 =", G0, ", alpha.0 =", A0, "|", NSEED, "seeds\n\n")

res <- list(); k <- 0L
for (n in NS) {
  for (s in seq_len(NSEED)) {
    d  <- dat_sim(n = n, n.visit = 5, gamma.0 = G0, gamma.L = GL, alpha.0 = A0,
                  alpha.A = -0.04, alpha.L = 0.015, alpha.U = 0.015, seed = 5000 + s)
    ds <- ste_prepare(d, n.visit = 5, truncate = 0.95)
    w_real <- ds$ipw.stab.cum
    # the control, built on the SAME L values the propensity actually uses
    w_nc <- 1 + exp(G0 + GL * ds$L.lead1[is.finite(ds$L.lead1)])
    for (nm in c("real", "control")) {
      w <- if (nm == "real") w_real else w_nc
      w <- w[is.finite(w) & w > 0]
      k <- k + 1L
      res[[k]] <- data.frame(
        n = n, seed = s, which = nm, nw = length(w),
        hill_1pct = hill(w, max(10L, floor(0.01 * length(w)))),
        hill_2pct = hill(w, max(10L, floor(0.02 * length(w)))),
        hill_5pct = hill(w, max(10L, floor(0.05 * length(w)))),
        ess_frac  = ess_frac(w), max_w = max(w),
        me_lo = mean_excess(w, quantile(w, 0.90)), me_hi = mean_excess(w, quantile(w, 0.99)))
    }
  }
  cat("  n =", n, "done\n")
}
r <- do.call(rbind, res)
saveRDS(r, file.path(ROOT, "dev/negative_control_results.rds"))

cat("\n1. HILL TAIL-INDEX ESTIMATES (median over seeds).\n")
cat("   alpha_hat <= 2 is conventionally read as 'infinite variance'.\n")
cat("   THE CONTROL HAS alpha = INFINITY BY ALGEBRA -- any finite estimate for it is a FALSE POSITIVE.\n")
print(as.data.frame(r %>% group_by(n, which) %>%
  summarise(hill_1pct = round(median(hill_1pct, na.rm = TRUE), 2),
            hill_2pct = round(median(hill_2pct, na.rm = TRUE), 2),
            hill_5pct = round(median(hill_5pct, na.rm = TRUE), 2), .groups = "drop")))

cat("\n2. ESS/n AND MAX-WEIGHT SCALING IN n (median over seeds).\n")
sc <- r %>% group_by(n, which) %>%
  summarise(ess_frac = median(ess_frac), max_w = median(max_w), .groups = "drop")
print(as.data.frame(sc %>% mutate(across(where(is.numeric), ~signif(., 4)))))
cat("\n   scaling exponents (slope of log(y) on log(n)):\n")
for (nm in c("real", "control")) {
  x <- sc %>% filter(which == nm)
  cat(sprintf("     %-8s  d log(ess_frac)/d log n = %+.3f   d log(max_w)/d log n = %+.3f\n", nm,
              coef(lm(log(x$ess_frac) ~ log(x$n)))[2], coef(lm(log(x$max_w) ~ log(x$n)))[2]))
}

cat("\n3. MEAN EXCESS (median over seeds) -- rising e(u) is read as a heavy tail.\n")
print(as.data.frame(r %>% group_by(n, which) %>%
  summarise(e_at_q90 = round(median(me_lo, na.rm = TRUE), 4),
            e_at_q99 = round(median(me_hi, na.rm = TRUE), 4),
            rising = ifelse(median(me_hi, na.rm = TRUE) > median(me_lo, na.rm = TRUE), "yes", "no"),
            .groups = "drop")))

cat("\n#### VERDICT ####\n")
cat("If the CONTROL rows show finite Hill estimates, a collapsing ESS/n, a growing max, and rising\n")
cat("mean excess -- i.e. the same qualitative signature as the real weights -- then all four\n")
cat("diagnostics FIRE ON A VARIABLE WITH NO TAIL INDEX AT ALL. They therefore cannot establish a\n")
cat("tail-index/infinite-variance claim for the real weights, and the tail-index mechanism is dead\n")
cat("on evidence rather than on assertion.\n")
cat("\nNEGATIVE CONTROL DONE (saved dev/negative_control_results.rds)\n")
