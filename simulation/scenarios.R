# Scenario grid for the ADEMP simulation study.
# Four axes x 3 levels = 81 cells, mirroring Limozin, Seaman & Su (2025) but
# realised through Keogh's additive-hazard DGM (steCI::dat_sim):
#   - sample size n
#   - event rate       via alpha.0 (low/med/high); low probes the negative-
#                      hazard regime that stresses additive hazards
#   - confounding      via gamma.L (L -> A strength)
#   - treatment prev.  via gamma.0 (intercept of the treatment model)
# True marginal survival / risk difference depends only on alpha.0 (the
# treatment regime is fixed in the pseudo-RCT), so there are 3 distinct truths.

make_scenarios <- function() {
  # Event-rate levels are constrained by hazard non-negativity: an additive
  # hazard cannot be pushed very low without h = alpha.0 + alpha.L*L + ... going
  # negative for the treated arm (whose L drifts down). alpha.0 = 0.05 gave a
  # ~44% negative-hazard truncation and a degenerate DGM; 0.15/0.25/0.40 keep
  # the truncation negligible (threshold L < -(alpha.0-0.04)/0.015). See
  # simulation/calibrate_dgm.R for the calibration.
  g <- expand.grid(
    n       = c(300L, 1000L, 5000L),
    alpha.0 = c(0.15, 0.25, 0.4),  # event rate: low / med / high (achievable range)
    gamma.L = c(0.5, 1.5, 3.0),    # confounding: weak / moderate / strong
    gamma.0 = c(0.0, -1.0, -2.0),  # prevalence: high / med / low
    stringsAsFactors = FALSE)
  g$event_rate <- c("low", "med", "high")[match(g$alpha.0, c(0.15, 0.25, 0.4))]
  g$prevalence <- c("high", "med", "low")[match(g$gamma.0, c(0.0, -1.0, -2.0))]
  g$confounding <- c("weak", "moderate", "strong")[match(g$gamma.L, c(0.5, 1.5, 3.0))]
  g$scenario <- sprintf("n%d_ev%s_cf%s_pv%s", g$n, g$event_rate,
                        g$confounding, g$prevalence)
  g$scen_id <- seq_len(nrow(g))
  g
}

# Fixed DGM parameters (shared across all scenarios).
DGM_FIXED <- list(n.visit = 5, alpha.A = -0.04, alpha.L = 0.015, alpha.U = 0.015)

# Analysis settings.
T_HOR <- 1:5
