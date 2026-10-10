# The TRUE Jensen split on the main run.
# CORRECTION: se_ratio_fold = mean_se/fold_sd changes the DENOMINATOR (tail-robust scale).
# The JENSEN-free measure is rms_se/emp_sd, which changes the NUMERATOR:
#   E[se^2] = mean_se^2 (1+CV^2)  =>  se_ratio = (rms_se/emp_sd) * 1/sqrt(1+CV^2).
# An unbiased variance estimator has E[se^2] = Var(est), i.e. rms_se/emp_sd = 1. These are different
# corrections and must not be conflated.
suppressWarnings(suppressMessages({library(dplyr); library(tidyr)}))
ROOT <- "."
NEW  <- file.path(ROOT, "data/results")
source(file.path(ROOT, "simulation/scenarios.R"))
scn <- make_scenarios()[, c("scenario","n","gamma.L","gamma.0","alpha.0")]
rp  <- readRDS(file.path(NEW, "replicates_all.rds"))

pt <- rp %>% filter(method == "point", estimand == "mrd", is.finite(estimate), abs(estimate) <= 1) %>%
  group_by(scenario, tau) %>% summarise(emp_sd = sd(estimate), .groups = "drop")
se <- rp %>% filter(method == "m6_analytic_cluster", estimand == "mrd",
                    is.finite(se), se > 0, se <= 0.5) %>%
  group_by(scenario, tau) %>%
  summarise(mean_se = mean(se), rms_se = sqrt(mean(se^2)), cv_se = sd(se)/mean(se), .groups = "drop")
d <- pt %>% inner_join(se, by = c("scenario","tau")) %>% inner_join(scn, by = "scenario") %>%
  mutate(published = mean_se/emp_sd, jensen_free = rms_se/emp_sd,
         jensen_factor = 1/sqrt(1+cv_se^2))

cat("######## TRUE JENSEN SPLIT -- MAIN RUN ########\n\n")
cat("identity check |published - jensen_free*jensen_factor| max =",
    format(max(abs(d$published - d$jensen_free*d$jensen_factor), na.rm = TRUE), digits = 3), "\n\n")

cat("At gamma.L=3, tau=5 (the headline cell), by n:\n")
h <- d %>% filter(gamma.L == 3, tau == 5) %>% group_by(n) %>%
  summarise(cv_se = round(mean(cv_se),3), published = round(mean(published),3),
            jensen_free = round(mean(jensen_free),3),
            jensen_share = round(mean(-log(jensen_factor))/mean(-log(published)),3), .groups="drop")
print(as.data.frame(h))

cat("\nAveraged over ALL 81 scenarios, by tau (n=5000):\n")
a <- d %>% filter(n == 5000) %>% group_by(tau) %>%
  summarise(cv_se = round(mean(cv_se),3), published = round(mean(published),3),
            jensen_free = round(mean(jensen_free),3), .groups="drop")
print(as.data.frame(a))

cat("\nThe n-trend, Jensen-free, at gamma.L=3 tau=5:\n")
t <- h %>% select(n, published, jensen_free)
cat("  published  :", paste(sprintf("%.3f", t$published), collapse=" -> "),
    sprintf(" (drop %.3f)", t$published[1]-t$published[nrow(t)]), "\n")
cat("  jensen-free:", paste(sprintf("%.3f", t$jensen_free), collapse=" -> "),
    sprintf(" (drop %.3f)", t$jensen_free[1]-t$jensen_free[nrow(t)]), "\n")
cat("\nJENSEN SPLIT DONE\n")
