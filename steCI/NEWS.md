# steCI 0.0.1

Released with the revised paper *Confidence Intervals and Simultaneous Bands for
Additive-Hazards Sequential Trial Emulation*.

## New features

* New `ties` argument (`"timereg"`, the default, or `"simultaneous"`) for
  `ste_fit()`, `ste_estimate()`, `ci_analytic()`, `ci_bootstrap()` and
  `ci_jackknife()`. `timereg::aalen()` breaks tied event times by random jitter,
  so a fit on tied data consumes the random-number stream. With
  `ties = "simultaneous"`, tied events enter one increment on a common risk set:
  the constant hazard difference and its model-based, robust and cluster-robust
  variances come from the new internal `.semiaalen()` (the same semiparametric
  model as timereg, including its `max.clust = 1000` cluster binning), and the
  cumulative coefficients for the risk difference come from `.aalen_if()`. On data
  without tied event times the two settings agree to machine precision.
* `.aalen_if()` gains `influence = FALSE` (cumulative coefficients only) and
  `singular = c("error", "skip")`. Under `"skip"` the contribution of a risk set
  whose weighted Gram matrix cannot be inverted is set to zero, and the count is
  returned as `n_singular`.
* `ci_bootstrap()` and `ci_jackknife()` accept `point =`, a full-sample
  `ste_estimate()` result to centre the intervals on. Its tie convention must
  match `ties`. Under `ties = "simultaneous"`, `ci_bootstrap()` also reports
  `n_singular`, the number of refits that skipped a singular risk set.
* New `time.grid` argument for `ste_prepare()`, `ste_estimate()`, `ci_bootstrap()`,
  `ci_jackknife()` and `ci_lef()`. For data recorded on a discrete sub-visit grid
  (e.g. `time.grid = 30` for daily data on 30-day visits), it snaps the
  trial-relative stop times to the exact grid, so events on the same
  trial-relative day are exact ties rather than last-bit near-ties.

## Unchanged

* The defaults (`ties = "timereg"`, `time.grid = NULL`) reproduce the published
  simulation replicates (every method) to within 1e-12
  (`dev/heart_ties/verify_default_path.R`).
* The applied Stanford heart-transplant example now uses
  `ties = "simultaneous"` and `time.grid = 30` throughout (`applied/heart_common.R`).

# steCI 0.0.0.9000

* Initial public release with the reproducibility compendium.
