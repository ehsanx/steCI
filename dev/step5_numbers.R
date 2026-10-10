# The exact numbers Step 5 needs, all from the canonical aggregation.
# Prints only what the rewritten sentences will quote, so nothing is typed from memory.
#   Rscript dev/step5_numbers.R
suppressWarnings(suppressMessages({library(dplyr); library(tidyr)}))
# Coverages depend on the truth, so the default is the exact-truth aggregation the paper
# reports (data/exact/results); RESULTS_DIR=data/results gives the earlier Monte-Carlo-truth figures.
ROOT <- Sys.getenv("REPO_ROOT", ".")
RES  <- Sys.getenv("RESULTS_DIR", file.path(ROOT, "data/exact/results"))
source(file.path(ROOT, "simulation/scenarios.R"))
ci  <- read.csv(file.path(RES, "metrics_ci.csv"), stringsAsFactors = FALSE)
scn <- make_scenarios()[, c("scenario", "n", "gamma.L", "gamma.0", "alpha.0")]
ci  <- merge(ci, scn, by = "scenario")
cat("SOURCE:", RES, "\n\n")

cat("A. PER-HORIZON coverage, mrd, mean over all 81 scenarios (for the L805 rewrite)\n")
a <- ci %>% filter(estimand == "mrd",
                   method %in% c("m6_analytic_cluster","m6_multiplier","boot_percentile","m6_band")) %>%
  group_by(method, tau) %>% summarise(cov = round(mean(coverage), 3), .groups = "drop") %>%
  pivot_wider(names_from = tau, values_from = cov, names_prefix = "tau")
print(as.data.frame(a))

cat("\nB. BIAS-ELIMINATED vs RAW coverage (is the residual bias or variance?), mrd, by n\n")
b <- ci %>% filter(estimand == "mrd", method == "m6_analytic_cluster") %>%
  group_by(n) %>% summarise(coverage = round(mean(coverage), 3),
                            coverage_be = round(mean(coverage_be), 3),
                            gain = round(mean(coverage_be) - mean(coverage), 4), .groups = "drop")
print(as.data.frame(b))
cat("   (E1 comparison -- where the bias story DOES hold)\n")
b1 <- ci %>% filter(estimand == "hd", method == "cluster") %>%
  group_by(n) %>% summarise(coverage = round(mean(coverage), 3),
                            coverage_be = round(mean(coverage_be), 3),
                            gain = round(mean(coverage_be) - mean(coverage), 4), .groups = "drop")
print(as.data.frame(b1))

cat("\nC. SIMULTANEOUS (joint) coverage -- all methods (for the band cluster)\n")
j <- read.csv(file.path(RES, "metrics_joint.csv"), stringsAsFactors = FALSE)
jj <- j %>% filter(estimand == "mrd") %>% group_by(method) %>%
  summarise(joint = round(weighted.mean(joint_coverage, n_curves), 4),
            curves = sum(n_curves), .groups = "drop") %>% arrange(desc(joint))
print(as.data.frame(jj))

cat("\nD. BAND: pointwise vs simultaneous, by n (the abstract/L215 claim)\n")
bp <- ci %>% filter(estimand == "mrd", method == "m6_band") %>%
  group_by(n) %>% summarise(pointwise = round(mean(coverage), 3), .groups = "drop")
bj <- j %>% filter(estimand == "mrd", method == "m6_band") %>%
  merge(scn, by = "scenario") %>% group_by(n) %>%
  summarise(simultaneous = round(weighted.mean(joint_coverage, n_curves), 3),
            curves = sum(n_curves), .groups = "drop")
print(as.data.frame(merge(bp, bj, by = "n")))

cat("\nE. Band pointwise by horizon (is 'pointwise exceeds nominal' true at every tau?)\n")
e <- ci %>% filter(estimand == "mrd", method == "m6_band") %>%
  group_by(tau) %>% summarise(cov = round(mean(coverage), 3), .groups = "drop")
print(as.data.frame(e))
cat("\nSTEP5 NUMBERS DONE\n")
