# PHASE A (remediation plan): measure the truth from artefacts already on disk.
# NO re-run, NO ARC. Read-only except for writing dev/phaseA_results.rds.
#   Rscript dev/phaseA_measure.R
suppressWarnings(suppressMessages(library(dplyr)))
RES <- "./simulation/output/results"
ci  <- read.csv(file.path(RES, "metrics_ci.csv"), stringsAsFactors = FALSE)
tc  <- readRDS("./simulation/output/truth_cache.rds")
rp  <- readRDS(file.path(RES, "replicates_all.rds"))

# truth lookup mirroring the pipeline: keyed by alpha.0, indexed by tau
truth_mrd <- function(a0, tau) vapply(seq_along(a0), function(i)
  tc[[as.character(a0[i])]]$mrd[tau[i]], numeric(1))

cat("################ PHASE A RESULTS ################\n\n")

## ---- A1. SIMULTANEOUS band coverage (the headline claim, never measured) ----
b <- rp %>% filter(estimand == "mrd", method == "m6_band",
                   !is.na(lower), !is.na(upper), !is.na(tau))
b$truth <- truth_mrd(b$alpha.0, b$tau)
b$cov_pt <- b$lower <= b$truth & b$upper >= b$truth

# validation: pointwise average must reproduce the published ~0.960
cat("A1a. band POINTWISE coverage (validation vs published 0.961/0.959/0.961):\n")
print(b %>% group_by(n) %>% summarise(pointwise = mean(cov_pt), .groups = "drop"))
cat("    overall pointwise =", round(mean(b$cov_pt), 4), "\n\n")

# simultaneous: a curve counts only if all 5 horizons are present
sim <- b %>% group_by(scenario, n, rep) %>%
  summarise(k = n(), all_cov = all(cov_pt), .groups = "drop") %>% filter(k == 5)
cat("A1b. band SIMULTANEOUS coverage (all 5 horizons at once):\n")
cat("    curves =", nrow(sim), "| overall =", round(mean(sim$all_cov), 4),
    "| MC SE =", round(sqrt(mean(sim$all_cov) * (1 - mean(sim$all_cov)) / nrow(sim)), 4), "\n")
print(sim %>% group_by(n) %>% summarise(simultaneous = round(mean(all_cov), 4), .groups = "drop"))
by_scn <- sim %>% group_by(scenario) %>% summarise(s = mean(all_cov), .groups = "drop")
cat("    scenarios reaching 0.95:", sum(by_scn$s >= 0.95), "of", nrow(by_scn),
    "| min =", round(min(by_scn$s), 3), "| max =", round(max(by_scn$s), 3), "\n\n")

## ---- A2. per-horizon coverage (Table 4 averages this away) ----
cat("A2. mrd coverage by horizon tau (mean over 81 scenarios):\n")
ph <- ci %>% filter(estimand == "mrd",
                    method %in% c("m6_band","m6_multiplier","m6_analytic_cluster",
                                  "m6_analytic_row","boot_percentile")) %>%
  group_by(method, tau) %>% summarise(cov = round(mean(coverage), 3), .groups = "drop")
print(as.data.frame(tidyr::pivot_wider(ph, names_from = tau, values_from = cov,
                                       names_prefix = "tau")))
cat("\n")

## ---- A3. E2: is the residual bias or variance? ----
cat("A3. mrd raw vs bias-eliminated coverage and SE ratio, by n:\n")
scn <- rp %>% distinct(scenario, n)
e2 <- ci %>% filter(estimand == "mrd",
                    method %in% c("m6_analytic_cluster","m6_multiplier","m6_band","boot_percentile")) %>%
  left_join(scn, by = "scenario") %>%
  group_by(method, n) %>%
  summarise(cov = round(mean(coverage), 3), cov_be = round(mean(coverage_be), 3),
            gain_from_recentring = round(mean(coverage_be) - mean(coverage), 3),
            se_ratio = round(mean(se_ratio), 3), .groups = "drop")
print(as.data.frame(e2))
cat("\n    (E1 comparison: recentring gain for the cluster SE)\n")
e1 <- ci %>% filter(estimand == "hd", method == "cluster") %>% left_join(scn, by = "scenario") %>%
  group_by(n) %>% summarise(cov = round(mean(coverage), 3), cov_be = round(mean(coverage_be), 3),
                            gain = round(mean(coverage_be) - mean(coverage), 3),
                            se_ratio = round(mean(se_ratio), 3), .groups = "drop")
print(as.data.frame(e1))
cat("\n")

## ---- A4. failure taxonomy: lost compute vs genuine ----
cat("A4. replicate counts per scenario (18 cells are short):\n")
reps <- rp %>% group_by(scenario, scen_id, n, alpha.0, gamma.0, gamma.L) %>%
  summarise(n_rep = n_distinct(rep), .groups = "drop") %>% arrange(n_rep)
short <- reps %>% filter(n_rep < 1000)
cat("    short cells:", nrow(short), "of", nrow(reps), "\n")
short$prevalence <- c("high","med","low")[match(short$gamma.0, c(0, -1, -2))]
short$batch_multiple_of_32 <- short$n_rep %% 32 == 0
print(as.data.frame(short %>% select(scenario, n, prevalence, n_rep, batch_multiple_of_32)))
cat("\n    by n:\n"); print(as.data.frame(short %>% count(n, name = "n_short_cells")))
cat("    by prevalence:\n"); print(as.data.frame(short %>% count(prevalence, name = "n_short_cells")))

saveRDS(list(sim_band = sim, by_scn = by_scn, per_horizon = ph, e2 = e2, e1 = e1,
             reps = reps, short = short), "./dev/phaseA_results.rds")
cat("\nPHASE_A DONE (saved dev/phaseA_results.rds)\n")
