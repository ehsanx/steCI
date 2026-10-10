# Publication figures from the aggregated ADEMP metrics. Reads a metrics object
# (list with $ci, $empirical) from an .rds and writes PDFs to figures/.
# Defaults to the shipped full-run exact-truth metrics.
#   Rscript simulation/figures.R [path/to/metrics.rds]
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(timereg); library(survival)})
for (f in list.files("steCI/R", pattern = "[.]R$", full.names = TRUE)) source(f)
source("simulation/scenarios.R")

args <- commandArgs(trailingOnly = TRUE)
# Default is the CANONICAL exact metrics, not a pilot file that no longer exists.
# The previous default pointed at metrics_pilot_full.rds, which is absent -- so a bare
# `Rscript simulation/figures.R` failed rather than silently plotting the wrong data, but it
# also meant the default could not regenerate the shipped figures.
mpath <- if (length(args)) args[1] else
  "data/exact/results/metrics.rds"
if (!file.exists(mpath))
  stop("metrics file not found: ", mpath,
       "\n  pass one explicitly, or use the default data/exact/results/metrics.rds")
outdir <- Sys.getenv("FIGURES_OUT", "figures")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
m <- readRDS(mpath)
scn <- make_scenarios()[, c("scenario", "n", "event_rate", "prevalence", "confounding")]
ci <- merge(m$ci, scn, by = "scenario")

lab <- c(model = "Model (M1)", robust = "Robust, row (M2)", cluster = "Cluster-robust (M3)",
         boot_percentile = "Bootstrap (M4)", boot_basic = "Bootstrap-basic", boot_normal = "Bootstrap-normal",
         m6_analytic_row = "M6 analytic, row", m6_analytic_cluster = "M6 analytic, cluster",
         m6_multiplier = "M6 multiplier", m6_band = "M6 band (simult.)",
         m6a2_analytic_cluster = "M6-A2 cluster", m6a2_multiplier = "M6-A2 multiplier")
pal <- c("Model (M1)"="#999999","Robust, row (M2)"="#E69F00","Cluster-robust (M3)"="#0072B2",
         "Bootstrap (M4)"="#D55E00","M6 analytic, row"="#E69F00","M6 analytic, cluster"="#0072B2",
         "M6 multiplier"="#009E73","M6 band (simult.)"="#CC79A7","M6-A2 cluster"="#56B4E9",
         "M6-A2 multiplier"="#661100")
base_theme <- theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(),
                                               legend.title = element_blank())

## ---- Fig 1: E1 coverage & SE-ratio by method x n ----
hd_meths <- c("model", "robust", "cluster", "boot_percentile")
hd <- ci %>% filter(estimand == "hd", method %in% hd_meths) %>%
  group_by(n, method) %>% summarise(coverage = mean(coverage), se_ratio = mean(se_ratio), .groups = "drop") %>%
  mutate(Method = factor(lab[method], levels = lab[hd_meths]))
f1 <- hd %>% pivot_longer(c(coverage, se_ratio), names_to = "metric") %>%
  mutate(metric = recode(metric, coverage = "Coverage (nominal 0.95)", se_ratio = "SE ratio (target 1.0)")) %>%
  ggplot(aes(n, value, colour = Method, group = Method)) +
  geom_hline(data = data.frame(metric = c("Coverage (nominal 0.95)", "SE ratio (target 1.0)"), y = c(0.95, 1)),
             aes(yintercept = y), linetype = 2, colour = "grey40") +
  geom_line() + geom_point() + facet_wrap(~metric, scales = "free_y") +
  scale_x_log10(breaks = c(300, 1000, 5000)) + scale_colour_manual(values = pal) +
  labs(x = "sample size n", y = NULL,
       title = "E1: constant hazard difference — coverage and SE accuracy") + base_theme
ggsave(file.path(outdir, "fig1_e1_coverage.pdf"), f1, width = 9, height = 3.6)

## ---- Fig 2: E2 MRD coverage by method x horizon, one panel per n ----
mrd_meths <- intersect(c("m6_analytic_row","m6_analytic_cluster","m6_multiplier","m6_band",
                         "m6a2_analytic_cluster","m6a2_multiplier","boot_percentile"),
                       unique(ci$method[ci$estimand=="mrd"]))
mrd <- ci %>% filter(estimand == "mrd", method %in% mrd_meths) %>%
  group_by(tau, n, method) %>% summarise(coverage = mean(coverage), .groups = "drop") %>%
  mutate(Method = factor(lab[method], levels = lab[mrd_meths]))
f2 <- ggplot(mrd, aes(tau, coverage, colour = Method, group = Method)) +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "grey40") +
  geom_line() + geom_point() + facet_wrap(~n, labeller = label_both) +
  scale_colour_manual(values = pal) +
  labs(x = "horizon", y = "coverage",
       title = "E2: marginal risk difference — coverage by horizon and sample size") + base_theme
ggsave(file.path(outdir, "fig2_e2_coverage.pdf"), f2, width = 8, height = 4)

## ---- Fig 3: E1 coverage by event rate (where under-coverage bites) ----
f3 <- ci %>% filter(estimand == "hd", method %in% hd_meths) %>%
  group_by(event_rate, method) %>% summarise(coverage = mean(coverage), .groups = "drop") %>%
  mutate(Method = factor(lab[method], levels = lab[hd_meths]),
         event_rate = factor(event_rate, c("low","med","high"))) %>%
  ggplot(aes(event_rate, coverage, fill = Method)) +
  geom_hline(yintercept = 0.95, linetype = 2, colour = "grey40") +
  geom_col(position = position_dodge()) + coord_cartesian(ylim = c(0.6, 1)) +
  scale_fill_manual(values = pal) +
  labs(x = "event rate", y = "coverage", title = "E1 coverage by event rate") + base_theme
ggsave(file.path(outdir, "fig3_e1_by_eventrate.pdf"), f3, width = 7, height = 4)

## ---- Fig 4: example MRD curve with 95% simultaneous band ----
set.seed(20240704); d <- dat_sim(n = 1000, gamma.0 = -1, gamma.L = 1.5, alpha.0 = 0.25, seed = 704)
lf <- ci_lef(d, t.hor = 1:5, B = 500, approach = "a1", seed = 1)
tv <- readRDS("dev/truth_cache_exact.rds")[[as.character(0.25)]]
stopifnot(isTRUE(tv$exact), identical(tv$t.hor, 1:5))
mc <- lf$mrd[lf$mrd$type == "multiplier", ]; bd <- lf$band
df4 <- data.frame(tau = 1:5, mrd = mc$estimate, plo = mc$lower, phi = mc$upper,
                  blo = bd$lower, bhi = bd$upper, truth = tv$mrd)
f4 <- ggplot(df4, aes(tau, mrd)) +
  geom_ribbon(aes(ymin = blo, ymax = bhi, fill = "95% simultaneous band"), alpha = .25) +
  geom_errorbar(aes(ymin = plo, ymax = phi, colour = "95% pointwise (multiplier)"), width = .12) +
  geom_line(aes(colour = "estimate")) + geom_point(aes(colour = "estimate")) +
  geom_line(aes(y = truth, colour = "truth"), linetype = 2) +
  scale_colour_manual(values = c("estimate"="#0072B2","95% pointwise (multiplier)"="#0072B2","truth"="black")) +
  scale_fill_manual(values = c("95% simultaneous band"="#0072B2")) +
  labs(x = expression(horizon~tau), y = "marginal risk difference",
       title = "M6: risk-difference curve with pointwise CI and simultaneous band",
       subtitle = "n=1000, medium event rate/prevalence/confounding") + base_theme
ggsave(file.path(outdir, "fig4_mrd_band.pdf"), f4, width = 7, height = 4.4)

cat("Wrote:", paste(list.files(outdir, "[.]pdf$"), collapse = ", "), "to", outdir, "\n")
