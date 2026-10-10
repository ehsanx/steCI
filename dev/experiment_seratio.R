# Empirical SE-ratio limit and the (C5) within-subject cross-record share.
#
# Proposition 2 gives V_row(tau) = V_cl(tau) - Delta(tau) with Delta the aggregate
# within-subject cross-record term (positive under condition C5). Hence the row-to-cluster SE
# ratio has the limit
#     rho = SE_row / SE_cl = sqrt(V_row / V_cl) = sqrt(1 - Delta / V_cl),
# and Delta / V_cl = 1 - rho^2 is the fraction of the true (cluster) variance the
# row-level sandwich omits. This script estimates rho and Delta/V_cl for the E1
# constant hazard difference from the EXISTING full-design output (no re-run): each
# replicate already carries the row-level robust SE and the subject-clustered SE.
#   Rscript dev/experiment_seratio.R          # Web Table S6
# Needs the per-replicate record data/results/replicates_all.rds (91 MB, not in git; see the
# README's "Large files"). The SE ratios are truth-free, so data/exact/results/replicates_all.rds
# gives the same numbers. Env: SERATIO_REPLICATES (input), SERATIO_OUT (per-scenario .rds).
suppressMessages({library(dplyr); library(tidyr)})
REP <- Sys.getenv("SERATIO_REPLICATES", "data/results/replicates_all.rds")
if (!file.exists(REP)) stop("per-replicate record not found: ", REP, " (see the README, Large files)")
r <- readRDS(REP)
r <- r[r$estimand == "hd", ]
r$se <- as.numeric(r$se); r$estimate <- as.numeric(r$estimate)
r$n <- as.integer(r$n); r$rep <- as.integer(r$rep)

# wide: one row per (scenario, rep) with the point estimate and the three SEs
pt <- r %>% filter(method == "point") %>%
  select(scenario, scen_id, n, alpha.0, rep, gamma = estimate)
se <- r %>% filter(method %in% c("model", "robust", "cluster")) %>%
  select(scenario, rep, method, se) %>%
  pivot_wider(names_from = method, values_from = se)
w <- inner_join(pt, se, by = c("scenario", "rep")) %>%
  filter(is.finite(gamma), is.finite(robust), is.finite(cluster), cluster > 0) %>%
  mutate(Vrow = robust^2, Vcl = cluster^2,
         rho = sqrt(Vrow / Vcl), delta_frac = 1 - Vrow / Vcl)

# per-scenario summaries (mirror tab:e1: empirical SD of gamma-hat, mean SEs)
scn <- w %>% group_by(scen_id, scenario, n, alpha.0) %>%
  summarise(sd_gamma   = sd(gamma),
            se_row_bar = mean(robust), se_cl_bar = mean(cluster),
            rho_hat    = sqrt(mean(Vrow) / mean(Vcl)),
            delta_frac = 1 - mean(Vrow) / mean(Vcl),
            .groups = "drop") %>%
  mutate(seratio_row = se_row_bar / sd_gamma,
         seratio_cl  = se_cl_bar  / sd_gamma)

# per-n (average over the 27 scenarios), and by event rate
by_n <- scn %>% group_by(n) %>%
  summarise(rho_hat = mean(rho_hat), delta_frac = mean(delta_frac),
            seratio_row = mean(seratio_row), seratio_cl = mean(seratio_cl),
            rho_from_table = mean(seratio_row) / mean(seratio_cl), .groups = "drop")
by_ev <- scn %>% mutate(ev = c("low","med","high")[match(alpha.0, c(0.15,0.25,0.40))]) %>%
  group_by(ev) %>% summarise(rho_hat = mean(rho_hat), delta_frac = mean(delta_frac), .groups="drop")

cat("=== SE-ratio limit rho = sqrt(V_row/V_cl) and (C5) share Delta/V_cl, by n ===\n")
print(as.data.frame(by_n), digits = 3)
cat("\n=== by event rate (averaged over n and the other axes) ===\n")
print(as.data.frame(by_ev), digits = 3)
cat("\nOverall: rho_hat =", round(mean(scn$rho_hat), 3),
    " Delta/V_cl =", round(mean(scn$delta_frac), 3),
    " (min Delta/V_cl over scenarios =", round(min(scn$delta_frac), 3),
    "-> C5 Delta>0 holds in", sum(scn$delta_frac > 0), "/", nrow(scn), "scenarios)\n")
OUT <- Sys.getenv("SERATIO_OUT", file.path(tempdir(), "seratio_by_scenario.rds"))
saveRDS(scn, OUT); cat("wrote", OUT, "\n")
cat("\nEXPERIMENT_SERATIO DONE\n")
