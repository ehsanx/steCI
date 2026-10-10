# dev/ — development & validation scripts

Not part of the `steCI` package or the simulation study; these back specific
claims in the paper and are kept for reproducibility. The repository README maps each
Web Table to its script ("Paper results → scripts").

- **`validate_m6.R`** — validates the M6 machinery (`.aalen_if` + `ci_lef`):
  (a) on data without tied event times the custom weighted-Aalen fit reproduces
  `timereg::aalen` cumulative coefficients to ~machine precision, and (b) the M6
  cluster-robust analytic and multiplier SEs for the marginal risk difference agree with
  the nonparametric bootstrap SE. Backs the "M6 matches the bootstrap" claim.
- **`validate_a2.R`** — validates M6 Approach-2 (weight-estimation uncertainty):
  the analytic weight derivative `dw/dtheta` and the sensitivity `A(t)=dB/dtheta`
  match finite differences (~1e-10).
- **`approach2_paired.R`** — the size of the Approach-2 correction on COMMON VALID
  replicates (both standard errors pass the validity check), from the per-replicate
  record: Web Table S1 and the pooled figures (signed median -0.35%, range -3.4% to
  +4.4%, 36% of cells beyond 1% in absolute value, exactly zero at tau = 1). It first
  reproduces the earlier unpaired table of `a2_ratio_table.R` and the mean SEs of
  `metrics_ci.csv`. Needs `data/exact/results/replicates_all.rds`; stored run:
  `approach2_paired.log`.
- **`benchmark_timing.R`** — per-operation cost at n = 5000 (prepare, prepare + two
  fits, the M6 fit, `ci_lef` with B = 200), backing the computational-cost table; the
  M4 cost per replicate is extrapolated as B x (prepare + two fits). Stored runs:
  `benchmark_timing.log`.
- **`truncation_by_alpha0.R`** — the fraction of truncated (negative) hazard
  evaluations by baseline hazard: 14-43% across the nine prevalence x confounding cells
  at alpha_0 = 0.05, at most 0.2% at alpha_0 >= 0.15 (the latter from the main
  simulation's own per-replicate values). Stored run: `truncation_by_alpha0.log`.
- **`heart_grid_sensitivity.R`** — Web Table S9 (visit-grid ladder of the heart
  example) and the finest-rung jackknife/bootstrap corroboration. Uses the tie
  convention of the application (`simultaneous`, set in `applied/heart_common.R`);
  `HEART_TIES=timereg` overrides it. Stored output: `heart_grid_sensitivity_results.rds`.
- **`heart_ties/`** — `verify_default_path.R` reruns complete canonical simulation
  replicates (every method) through the current package with its defaults and compares
  them with the per-replicate record (to <= 1e-12); `run_canonical_checks.sh` runs it on
  three replicates together with `verify_canonical_repro.R`. Stored runs:
  `heart_ties/verify_default_path.log`.
- **`s5_rescore/`** — Web Table S5 (fixed vs random standardisation).
  `run_s5_reps.R` re-runs the replicate stream of `experiment_stdvar.R` (7 cells, 400
  replicates; 200 at n = 5000) and stores the per-replicate estimates in `reps/` (no truth
  inside); `score_s5.R` scores them against `dev/truth_cache_exact.rds` in seconds
  (`s5_summary.csv`, `s5_scored.rds`, `score_s5.log`: identical cluster-Wald coverage at
  every horizon; multiplier coverage differs by at most 0.0025 at a horizon and 0.0005
  averaged). `S5_SMOKE=1` reruns two replicates per cell.
- **`m4_ties/`** — Web Table S10, the exploratory comparison of bootstrap tie handling
  (see `m4_ties/README.md`).
- **`plim_closed_form.R`, `plim_seedcheck.R`, `plim_ablation.R`** — Web Appendix E: the
  population limit and origin-pooling component (Web Table S11; seven independent
  streams) and the attribution ablation (Web Table S12, all three arms scored against
  exact truths via `surv_exact(alpha.U = )` in `truth_closed_form.R`).

Run from the repo root, e.g. `Rscript dev/validate_m6.R`.
