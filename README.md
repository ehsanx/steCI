# steCI — Confidence Intervals for Sequential Trial Emulation with Additive Hazards

Reproducible research compendium and companion R package **`steCI`** for the methods
paper *Confidence Intervals and Simultaneous Bands for Additive-Hazards Sequential Trial
Emulation* (M. Ehsan Karim, 2026). An earlier version of the paper circulated as a preprint
under the title *Cluster on the Subject, Not the Record: ...*.

It derives, and evaluates the coverage of, confidence-interval procedures for the
**weighted-Aalen additive-hazards sequential-trial-emulation (STE)** estimator of a
causal marginal risk difference, addressing — for the additive-hazards outcome model —
the open problem stated by Limozin, Seaman & Su (2025): *"inference procedures for STE
based on additive hazard models warrant further research."*

> **Two estimands.** **E1** — a constant hazard difference; **E2** — the marginal
> risk-difference (MRD) curve over follow-up.

When the marginal hazard effect varies over follow-up, E1 is a design-weighted
constant summary that depends on the emulated-trial structure and analysis weights.

---

## Contents at a glance

| You want… | Go to |
|---|---|
| Install and use the estimator + CI methods | [`steCI/`](steCI/) · [Install](#install) · [Quick start](#quick-start) · [User manual](steCI/vignettes/steCI.Rmd) |
| The figures/tables and the code + data behind each | [Where each result comes from](#where-each-result-comes-from) |
| To re-run the whole simulation study | [Reproduce the study](#reproduce-the-study) |
| The applied illustration | [`applied/`](applied/) · [Reproduce the applied example](#reproduce-the-applied-example) |
| The exact-truth / sensitivity / robustness analyses | [`dev/`](dev/) · [the results map](#where-each-result-comes-from) |
| Which script reproduces a given Web Table or number | [Paper results → scripts](#paper-results--scripts) |
| What changed in this release | [`steCI/NEWS.md`](steCI/NEWS.md) |

---

## Background

- **Keogh et al. (2023, *Stat. Med.*)** give the additive-hazards STE estimator of the
  marginal risk difference and recommend the nonparametric bootstrap for confidence
  intervals, but **did not report confidence-interval coverage**.
- **Limozin et al. (2025, *SMMR*)** compare sandwich / bootstrap / LEF-bootstrap /
  jackknife intervals, but only for a **pooled-logistic** marginal structural model,
  explicitly leaving the additive-hazards case to future work.
- Because each subject is reused across several stacked emulated trials, the
  stacked records are dependent. The default row-level robust ("sandwich") standard
  error reported by `timereg::aalen` omits an aggregate within-subject cross-record
  term. Under the paper's positive sign condition (C5), this omission makes it
  **asymptotically anticonservative**. The term was positive throughout the assessed
  E1 settings and at E2 horizons 2–5, and empirically near zero at the first E2
  horizon. This package derives the estimator's influence functions. **Clustering
  the robust variance on the subject** recovers most of the E1 coverage deficit.
  For E2, the **multiplier / LEF method's cluster-robust analytic intervals** recover
  most of the row-level deficit, and its **multiplier intervals** recover about half.
  Coverage remains below nominal; the simultaneous band's measured coverage is
  **0.881**.

## Confidence-interval methods

| Method | Estimand | Basis | Status |
|--------|----------|-------|--------|
| M1 model-based SE       | E1                 | Aalen model variance                    | empirically anticonservative — avoid |
| M2 robust (row) SE      | E1                 | row-level sandwich (`timereg` default)  | inconsistent when the omitted term is nonzero; anticonservative under (C5) — avoid |
| M3 cluster-robust SE    | E1                 | sandwich clustered on subject           | **recommended (E1)** — recovers most of the deficit |
| M4 nonparametric bootstrap | E1, E2 (MRD)    | refit trials + weights + fit per resample | valid under standard conditions (not established here); full-refit comparator |
| M5 jackknife            | E1                 | delete-one-subject                      | valid — optional |
| M6 multiplier / LEF     | E2 (MRD)           | perturb pre-computed influence functions; pointwise + sup-*t* band | **recommended (E2)** — analytic intervals recover most of the deficit; multiplier intervals about half; band simultaneous coverage 0.881 |

## Install

```r
# install.packages("remotes")
remotes::install_github("ehsanx/steCI", subdir = "steCI", build_vignettes = TRUE)
vignette("steCI", package = "steCI")   # the user manual
```

The installable package lives in the [`steCI/`](steCI/) subfolder; the rest of the
repository reproduces the simulation study and the applied example. The user manual
([`steCI/vignettes/steCI.Rmd`](steCI/vignettes/steCI.Rmd)) walks through a complete
analysis of the public Stanford heart transplant data: trial emulation, both estimands,
the subject-clustered standard errors, and the multiplier bootstrap with its simultaneous
band. Every exported function also has a help page (`help(package = "steCI")`).

## Quick start

```r
library(steCI)
set.seed(1)

## simulate from the additive-hazard MSM (Keogh et al. 2021)
d <- dat_sim(n = 1000, alpha.0 = 0.25, gamma.0 = -1, gamma.L = 1.5)

## point estimates: constant hazard difference (E1) and the marginal
## risk-difference curve (E2), by empirical standardisation of a weighted Aalen fit
est <- ste_estimate(d, estimand = "both")

## confidence intervals
ci_analytic(ste_prepare(d))   # M1 model / M2 row-robust / M3 cluster-robust analytic SEs
ci_lef(d, approach = "a1")    # M6 multiplier/LEF: pointwise intervals + simultaneous band (E2)
ci_bootstrap(d, B = 200)      # M4 nonparametric bootstrap (reference)
ci_jackknife(d)               # M5 delete-one-subject jackknife
```

`ci_lef(..., approach = "a2")` additionally propagates weight-estimation uncertainty
through the two logistic weight models (a first-order correction whose asymptotic validity
the paper does not establish; on common valid replicates it changes the risk-difference
standard errors by a signed median of -0.35%, range -3.4% to +4.4%).
`ci_lef(..., standardization = "random")` adds the sampling variability of the
standardisation population (the unconditional variance); the default `"fixed"` conditions
on it.

### Tied event times

`timereg::aalen()` breaks tied event times by random jitter, so on tied data every fit
consumes random numbers. `ste_fit()`, `ste_estimate()`, `ci_analytic()`, `ci_bootstrap()` and
`ci_jackknife()` take `ties = "timereg"` (the default, used in the simulation study) or
`ties = "simultaneous"` (tied events enter one increment on a common risk set; used in the
applied example). On data without tied event times the two agree to machine precision.
For data recorded on a discrete sub-visit grid, `time.grid` (e.g. `30` for daily data on
30-day visits) makes events on the same trial-relative day exact ties:

```r
est <- ste_estimate(d_daily, n.visit = 6, t.hor = 1:6, ties = "simultaneous", time.grid = 30)
ci_bootstrap(d_daily, n.visit = 6, t.hor = 1:6, B = 1000, seed = 1,
             ties = "simultaneous", time.grid = 30, point = est)
```

Bootstrap resamples always contain ties (duplicated subjects share event times); see
`dev/m4_ties/` for how much the convention moves the M4 standard errors in the simulation.

## Repository layout

```
steCI/            the R package (estimator + all CI methods + dat_sim + tests + vignette)
  R/  man/  tests/  vignettes/  DESCRIPTION  NAMESPACE  LICENSE

simulation/       the ADEMP simulation study (81-cell factorial design)
  scenarios.R       the factorial grid (n x event rate x confounding x prevalence)
  run_sim.R  one_replicate.R    run one replicate / drive the design
  aggregate.R  metrics.R        aggregate chunks -> ADEMP metrics (bias, coverage, SE ratio, ...)
  figures.R                     regenerate every manuscript figure (writes to figures/)
  report.R  calibrate_dgm.R  smoke_test.R
  slurm/            UBC ARC (Sockeye) Slurm array-job scripts + submit_*.sh launchers

data/             the SMALL aggregated results the tables/figures are computed from
  exact/results/    AUTHORITATIVE aggregation, scored against the exact closed-form truths
  results/          companion aggregation (truth-free quantities such as SE ratios)

dev/              analysis, validation, sensitivity & robustness scripts backing specific
                  manuscript claims, plus their small result caches (*.rds) and the exact
                  closed-form truth cache (truth_cache_exact.rds). See the results map below.

applied/          Stanford heart-transplant illustration (public survival::jasa data)
figures/          the generated figures (fig1_e1_coverage.pdf, ...); regenerate via figures.R
```

The large raw simulation output and the local `simulation/output/` working directory are
**not** tracked; they are regenerated by the pipeline. The small aggregated metrics needed
to reproduce every figure and table *are* included, under `data/`. One exception: the
per-replicate output of the robustness studies (Web Tables S4 and S7; 357 chunk files plus
the scored summary, about 0.9 MB) is shipped in `simulation/output/robust/`, with its Slurm
logs in `simulation/output/logs/robust_*`, so those tables can be re-scored without a
cluster run.

### Large files

The larger research inputs are not tracked in git, to keep the code repository small; they
are attached to GitHub release v0.0.1.
Download [steCI-research-inputs-0.0.1.zip](https://github.com/ehsanx/steCI/releases/download/v0.0.1/steCI-research-inputs-0.0.1.zip)
from that release and extract it at the repository root. It supplies the
paths below and a `SHA256SUMS.txt` manifest. These inputs are needed only by the listed
scripts; using the estimator, the manual and the applied example does not require them.

| File | Size | Needed by |
|---|---|---|
| `data/exact/results/replicates_all.rds` and the companion `data/results/replicates_all.rds` | 91 MB and 93 MB | `dev/approach2_paired.R`, `dev/truncation_by_alpha0.R` (part 2), `dev/experiment_seratio.R`, `dev/heart_ties/verify_default_path.R`, `dev/verify_canonical_repro.R`, `dev/check_maxclust.R`, `dev/m4_ties/extract_canonical.R` |
| `dev/m4_ties/reps/{pilot,prod}/` | 837 files, 12 MB | `dev/m4_ties/se_subset.R` (Web Table S10) |
| `dev/m4_ties/canon_boot.rds` | 32 MB | `dev/m4_ties/se_subset.R`; rebuilt from `replicates_all.rds` by `dev/m4_ties/extract_canonical.R` |

Every script that reads these inputs states the path at its top. Stored summaries
(`*.log`, `*.rds`, `*.csv`) are tracked where available. Re-scoring Web Table S6 and
re-running the paired Approach-2 and tie-handling checks require the companion archive.

### The two aggregations in `data/`

Both are the **same replicates**, re-aggregated two ways:

- **`data/exact/results/`** — scored against the **exact closed-form marginal truths**.
  This is authoritative for anything that depends on a truth (bias, MSE, coverage).
  `simulation/figures.R` reads `data/exact/results/metrics.rds` by default.
- **`data/results/`** — the companion aggregation; use it for **truth-free** quantities
  (e.g. empirical-vs-model SE ratios), which are identical across the two.
- **`data/truth_cache.rds`** — the single Monte-Carlo truth draw the main run was first scored
  against (superseded by `dev/truth_cache_exact.rds`; kept because `dev/truth_closed_form.R`,
  `dev/verify_closed_form.R` and the truth-rebuild scripts compare against it).

## Where each result comes from

Each row: the reported object → the script that produces it → the input it reads.
(Paths are relative to the repo root; run every script from there. The scripts listed under
[Paper results → scripts](#paper-results--scripts) default to the shipped locations, and their
environment overrides (`REPO_ROOT`, `RESULTS_DIR` or script-specific ones) are documented at the
top of each script. Some older exploratory `dev/` scripts (`phaseA_*`, `phaseC_*`, the ARC
collectors) read the local working directory `simulation/output/`, which exists only after
re-running the pipeline; the exception is `dev/robust_collect.R`, whose input
`simulation/output/robust/` is shipped.)

| Result | Script | Reads |
|---|---|---|
| **Figures** (E1 coverage, E2 coverage, MRD band, …) | `simulation/figures.R` | `data/exact/results/metrics.rds` → writes `figures/` |
| **Main ADEMP metrics** (bias, coverage, SE ratio by method × cell) | `simulation/aggregate.R` + `simulation/metrics.R` | regenerates `data/…` from per-replicate chunks |
| **Exact closed-form truths** (E1 hazard diff., E2 MRD curve) | `dev/truth_closed_form.R` | → `dev/truth_cache_exact.rds`; `dev/rebuild_truth_multiseed.R` is a historical Monte Carlo comparison |
| **`max.clust` sensitivity** (M3 cluster SE under exact vs binned clustering) | `dev/maxclust_arc.R` (runner) → `dev/maxclust_collect.R` (scoring) | → `dev/maxclust_results.rds` |
| **Population-limit / finite-n bias of the Aalen fit** | `dev/plim_arc.R`, `dev/plim_closed_form.R` → `dev/plim_collect.R` | `data/exact/results/`, → `dev/plim_results.rds` |
| **Jensen decomposition of the coverage shortfall** (exploratory; backs no paper result) | `dev/jensen_decomposition.R` | `data/results/`, → `dev/jensen_decomposition_results.rds` |
| **Weight-truncation sensitivity** (exploratory; backs no paper result; unrelated to the hazard-truncation fractions below) | `dev/truncation_arc.R` (ARC runner, `simulation/slurm/truncation.slurm`) → `dev/truncation_collect.R`, `dev/truncation_paired.R`, `dev/truncation_allcaps.R`; local runs `dev/truncation_redo.R`, `dev/probe_truncation_power.R` | `simulation/output/truncation/` (not shipped); `truncation_collect.R`, `truncation_paired.R`, `truncation_allcaps.R` and `probe_truncation_power.R` score against a superseded Monte-Carlo truth cache, `simulation/output/truth_cache.rds` (not shipped) → `dev/truncation_summary.rds`, `dev/truncation_paired_results.rds`, `dev/truncation_allcaps.rds` |
| **Negative-control diagnostic** (exploratory; backs no paper result) | `dev/negative_control.R` | → `dev/negative_control_results.rds` |
| **Robustness (misspecification, time-varying effect, 2nd confounder)** | `dev/robustness_*.R`, `dev/robust_arc.R` → `dev/robust_collect.R` | `simulation/output/robust/` (shipped), `dev/robustness_*.rds` → `simulation/output/robust/robust_scored.rds` |
| **Variance-formula / M6 / Approach-2 validation** | `dev/validate_*.R`, `dev/verify_*.R` | (self-contained checks) |
| **Applied Stanford heart-transplant analysis** | `applied/heart_analysis.R` (+ `applied/heart_common.R`) | `survival::jasa` (public) → `applied/heart_results.rds` |
| **Approach-2 SE change on common valid replicates** | `dev/approach2_paired.R` | `data/exact/results/` (incl. `replicates_all.rds`) → `dev/approach2_paired.log` |
| **Bootstrap tie handling (M4)** | `dev/m4_ties/` (see its README) | canonical replicates + `dev/m4_ties/reps/` → `dev/m4_ties/se_subset.*` |
| **Computational cost** | `dev/benchmark_timing.R` | `dat_sim()` → `dev/benchmark_timing.log` |
| **Hazard-truncation fractions by event rate** | `dev/truncation_by_alpha0.R` | `dat_sim()` (+ `data/results/replicates_all.rds`) → `dev/truncation_by_alpha0.log` |

Paper-backing `dev/` scripts are named by the claim they back (e.g.
`dev/failure_rates_and_nominal.R`, `dev/truncation_by_alpha0.R`); where stored, their `*.rds`
caches or `*.log` runs hold the computed numbers.

## Paper results → scripts

Run from the repository root (`Rscript <script>`; on Windows run Rscript from PowerShell or
cmd). Web Tables are numbered as in the paper's Supplementary Material (Online Resource 1).

| Paper result | Script | Stored output | Notes |
|---|---|---|---|
| Tables 3–5, Figures 1–3 and Web Figure S1 (coverage and SE ratios; Web Figure S1 is E1 coverage by event rate) | `simulation/figures.R`, `dev/check_tables_vs_metrics.R` | `figures/`, `data/exact/results/`, `dev/check_tables_vs_metrics.log` | the check recomputes every printed cell of Tables 3–5 from `data/exact/results/` (0 mismatches); Figure 2 is coverage by horizon, one panel per n |
| Section 4 per-horizon E2 coverage, bias-eliminated coverage (Cov.†), band simultaneous coverage 0.881 | `dev/step5_numbers.R` | — (prints in seconds) | reads `data/exact/results/` |
| Band simultaneous coverage across scenarios (0.754–0.956, 3 of 81 at 0.95); per-method construction-failure rates | `dev/failure_rates_and_nominal.R` | — (prints in seconds) | reads `data/exact/results/` |
| Web Table S2 (scenario grid) | `simulation/scenarios.R` | — | |
| Web Tables S3 and S8 (E2 coverage and SE ratio by horizon) | `dev/check_tables_vs_metrics.R` | `dev/check_tables_vs_metrics.log` | |
| Web Tables S4 and S7 (robustness: second confounder; misspecified and time-varying-effect mechanisms) | `dev/robustness_2confounder.R`, `dev/robustness_misspec.R`, `dev/robustness_tvfx.R`, run on ARC by `dev/robust_arc.R` (`simulation/slurm/robust.slurm`) → `dev/robust_collect.R`; exact targets `dev/robustness_truth_exact.R` (S4), `dev/robustness_truth_cox.R` and `dev/robustness_truth_tvfx.R` (S7) | `simulation/output/robust/` (per-replicate ARC output, 357 chunks, and `robust_scored.rds`), `simulation/output/logs/robust_*`, `dev/robustness_frozen_truths.rds`, `dev/robustness_truth_exact.rds`, `dev/truth_misspec_exact.rds`, `dev/truth_tvfx_exact.rds` | `Rscript dev/robust_collect.R` re-scores the shipped chunks in seconds (set `ROBUST_OUT` to a copy to leave the shipped `robust_scored.rds` untouched) |
| Web Table S6 (condition C5, row-to-cluster SE ratio) | `dev/experiment_seratio.R` | — | needs `data/results/replicates_all.rds` |
| Influence-function validations: `.aalen_if()` vs `timereg` and M6 vs bootstrap; Approach-2 derivatives vs finite differences; unconditional vs fixed standardisation on one dataset | `dev/validate_m6.R`, `dev/validate_a2.R`, `dev/validate_stdvar.R` | — | self-contained, minutes each |
| Table 6 (per-replicate cost): about 21 s, 15 s, 70 min | `dev/benchmark_timing.R [n]` | `dev/benchmark_timing.log` (3 runs, n = 5000) | the M4 time is extrapolated from the measured prepare-plus-fit |
| Event-rate constraint: 14–43% truncated at α0 = 0.05, at most 0.2% at α0 ≥ 0.15 | `dev/truncation_by_alpha0.R` | `dev/truncation_by_alpha0.log` | part 2 needs `data/results/replicates_all.rds` |
| Tables 7–8 and Section 5 (Application), incl. the five bootstrap refits with a singular risk set | `applied/heart_analysis.R` | `applied/heart_results.rds` | `ties = "simultaneous"`, `time.grid = 30` (`applied/heart_common.R`); about 15 min |
| Web Table S1 (Approach-2 SE change): median −0.35%, range −3.4% to +4.4% | `dev/approach2_paired.R` | `dev/approach2_paired.log` | common valid replicates; `dev/a2_ratio_table.R` is the earlier unpaired version |
| Web Table S5 (fixed vs random standardisation); multiplier differences at most 0.0025 per horizon, 0.0005 averaged | `dev/s5_rescore/score_s5.R` (scores the stored per-replicate estimates of `dev/s5_rescore/run_s5_reps.R`); `dev/experiment_stdvar.R` is the original one-pass version, `dev/validate_stdvar.R` the single-dataset check | `dev/s5_rescore/{s5_summary.csv,s5_scored.rds,score_s5.log}`, `dev/s5_rescore/reps/` | scored against `dev/truth_cache_exact.rds`; uses `ci_lef(standardization = )`'s random-standardisation term |
| Web Table S9 (visit-grid sensitivity of the heart example) | `dev/heart_grid_sensitivity.R` | `dev/heart_grid_sensitivity_results.rds` | same tie convention as the application; 8 workers |
| Web Table S10 (bootstrap tie handling) | `dev/m4_ties/se_subset.R` | `dev/m4_ties/se_subset.{rds,csv,log}`, `se_subset_table.tex` | needs the large M4 inputs; see `dev/m4_ties/README.md` |
| Web Table S11 (population limit, origin pooling, seven streams) | `dev/plim_closed_form.R`, `dev/plim_seedcheck.R` | `dev/plim_closed_form*.rds`, `dev/plim_seedcheck.rds` | |
| Web Table S12 (attribution ablation) | `dev/plim_ablation.R` | `dev/plim_ablation.rds` | all three arms scored against exact truths (`surv_exact(alpha.U = )` in `dev/truth_closed_form.R`) |
| Exact closed-form truths | `dev/truth_closed_form.R`, `dev/verify_closed_form.R` | `dev/truth_cache_exact.rds` | |
| Package defaults unchanged by the tie switch | `dev/heart_ties/verify_default_path.R`, `dev/verify_canonical_repro.R` | `dev/heart_ties/verify_default_path.log` | needs `data/exact/results/replicates_all.rds` |

Exploratory only, backing no paper result: `dev/experiment_bca.R`, a probe of BCa and
related bootstrap intervals for the risk difference (two cells at n = 300). Its earlier
output was scored against a superseded single Monte-Carlo truth draw and was not stored;
the script now reads `dev/truth_cache_exact.rds` but has not been re-run. The paper does
not rely on it: it states only that interval shape has not been formally excluded as a
source of the residual shortfall.

## Reproduce the study

The shipped results can be checked without rerunning the simulation:

```sh
Rscript dev/check_tables_vs_metrics.R       # Tables 3–5: compare every printed cell
Rscript dev/step5_numbers.R                 # per-horizon and simultaneous coverage
Rscript dev/s5_rescore/score_s5.R            # re-score the shipped standardisation study
Rscript simulation/figures.R               # plot the shipped exact-truth metrics
Rscript simulation/smoke_test.R            # exercise the methods on small new datasets
```

The smoke test does not produce simulation chunks for aggregation. Exact closed-form
truths are included in `dev/truth_cache_exact.rds`; rebuild that cache with
`Rscript dev/truth_closed_form.R`.

To rerun the 81-cell, 1000-replicate design on a configured Slurm cluster:

```sh
bash simulation/slurm/submit_all.sh full
# The dependent aggregation writes simulation/output/results/ using exact truths.
# After it finishes, plot that new aggregation explicitly:
Rscript simulation/figures.R simulation/output/results/metrics.rds
```

See [`simulation/slurm/submit_all.sh`](simulation/slurm/submit_all.sh). Set your own
`--account`, `--mail-user` and module/library paths in the Slurm headers. Start each
profile with a fresh `simulation/output/intermediate/` directory; aggregation can
summarise partial runs, so check the completed replicate counts before interpreting
its output. The shipped canonical metrics in `data/exact/results/` are kept separately.

To aggregate existing chunks locally, use a separate R session:

```r
Sys.setenv(SIM_OUT = "simulation/output/intermediate",
           SIM_RESULTS = "simulation/output/results",
           SIM_TRUTH_CACHE = "dev/truth_cache_exact.rds")
source("simulation/aggregate.R")
```

`simulation/figures.R` writes to `figures/` by default. Set `FIGURES_OUT` to a different
directory when comparing regenerated figures with the shipped PDFs.

## Reproduce the applied example

```sh
Rscript applied/heart_analysis.R             # Stanford heart transplant; all six CI methods (~15 min)
Rscript dev/heart_grid_sensitivity.R         # Web Table S9 (visit-grid ladder; 8 PSOCK workers)
```

Both use `ties = "simultaneous"` and `time.grid = 30` for every printed heart number
(`applied/heart_common.R`); set the environment variable `HEART_TIES=timereg` to rerun them
under timereg's jittered ties. The seed is set before the first random draw, so reruns are
bit-identical within a setting.

## Data-generating mechanism (fixed)

Ported from Keogh et al. (2021, *Biometrical Journal*), with zero-based visits
`k = 0, ..., K-1`: an unmeasured frailty `U ~ N(0, 0.1^2)`;
`L_0 ~ N(U, 1)` and, for `k >= 1`, a time-varying confounder
`L_k ~ N(0.8 L_{k-1} - A_{k-1} + 0.1k + U, 1)`;
monotone treatment `logit Pr(A_k = 1 | A_{k-1}=0, L_k) = gamma_0 + gamma_L L_k`; and a
piecewise-constant generating hazard `h = alpha_0 + alpha_A A + alpha_L L + alpha_U U`
with `alpha_A = -0.04`, `alpha_L = alpha_U = 0.015`. Nonpositive values produce no
event in that interval, so the effective event hazard is `max(h, 0)`. The 81-cell grid varies `n`,
`alpha_0` (event rate), `gamma_L` (confounding), and `gamma_0` (prevalence).

## Citation

> Karim, M. E. *Confidence Intervals and Simultaneous Bands for Additive-Hazards
> Sequential Trial Emulation.* (2026.)

## References

- Keogh RH, Gran JM, Seaman SR, Davies G, Vansteelandt S (2023). Causal inference in
  survival analysis using longitudinal observational data: Sequential trials and
  marginal structural models. *Statistics in Medicine* 42(13):2191–2225.
- Limozin JM, Seaman SR, Su L (2025). Inference procedures in sequential trial
  emulation with survival outcomes. *Statistical Methods in Medical Research*.
- Keogh RH, Seaman SR, Gran JM, Vansteelandt S (2021). Simulating longitudinal data
  from marginal structural models using the additive hazard model. *Biometrical
  Journal* 63(7):1526–1541.
- Lin DY, Ying Z (1994). Semiparametric analysis of the additive risk model.
  *Biometrika* 81(1):61–71.

## License

MIT © 2026 M. Ehsan Karim. See [LICENSE](LICENSE).
