# Build the Approach-2 / Approach-1 SE-ratio table.
#
# WHY. The paper asserts the size of the weight-estimation correction in three places and points each
# time at Table 4 (tab:e2) -- but tab:e2 reports COVERAGE, not standard errors, so the claim has never
# actually been tabulated. This produces the numbers for a Web Table so the pointers can be redirected
# to something that contains the quantity.
#
# Canonical results only.
#   Rscript dev/a2_ratio_table.R
suppressWarnings(suppressMessages({library(dplyr); library(tidyr)}))
ROOT <- "."
RES  <- file.path(ROOT, "data/results")
source(file.path(ROOT, "simulation/scenarios.R"))
ci  <- read.csv(file.path(RES, "metrics_ci.csv"), stringsAsFactors = FALSE)
scn <- make_scenarios()[, c("scenario", "n")]
ci  <- merge(ci, scn, by = "scenario")

pair <- function(m1, m2, label) {
  a <- ci %>% filter(estimand == "mrd", method == m1) %>% select(scenario, n, tau, se1 = mean_se)
  b <- ci %>% filter(estimand == "mrd", method == m2) %>% select(scenario, tau, se2 = mean_se)
  inner_join(a, b, by = c("scenario", "tau")) %>%
    filter(is.finite(se1), is.finite(se2), se1 > 0) %>%
    mutate(pct = 100 * (se2 - se1) / se1, pair = label)
}
d <- bind_rows(pair("m6_analytic_cluster", "m6a2_analytic_cluster", "analytic cluster"),
               pair("m6_multiplier",       "m6a2_multiplier",       "multiplier"))

cat("######## APPROACH-2 vs APPROACH-1: % CHANGE IN THE MEAN SE ########\n")
cat("(canonical results; negative = A2 smaller)\n\n")

cat("BY HORIZON x SAMPLE SIZE, analytic cluster (median % change) -- the table body\n")
tb <- d %>% filter(pair == "analytic cluster") %>% group_by(n, tau) %>%
  summarise(med = round(median(pct), 2), .groups = "drop") %>%
  pivot_wider(names_from = tau, values_from = med, names_prefix = "tau")
print(as.data.frame(tb))

cat("\nBY HORIZON, pooled over n (median / max |%| / share exceeding 1%)\n")
th <- d %>% filter(pair == "analytic cluster") %>% group_by(tau) %>%
  summarise(median_pct = round(median(pct), 2), max_abs = round(max(abs(pct)), 2),
            frac_gt1 = round(mean(abs(pct) > 1), 3), .groups = "drop")
print(as.data.frame(th))

cat("\nOVERALL, per method pair\n")
ov <- d %>% group_by(pair) %>%
  summarise(cells = n(), median_pct = round(median(pct), 2),
            min_pct = round(min(pct), 2), max_pct = round(max(pct), 2),
            frac_gt1 = round(mean(abs(pct) > 1), 3), frac_gt2 = round(mean(abs(pct) > 2), 3),
            .groups = "drop")
print(as.data.frame(ov))

cat("\nSANITY: tau=1 should be ~0 -- the cumulative weights are identically 1 at the first visit,\n")
cat("so there is no weight-estimation uncertainty to propagate there.\n")
print(as.data.frame(d %>% filter(tau == 1) %>% group_by(pair) %>%
  summarise(max_abs_pct = signif(max(abs(pct)), 3), .groups = "drop")))
cat("\nA2 TABLE DONE\n")
