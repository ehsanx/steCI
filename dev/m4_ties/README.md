# dev/m4_ties/ — tie handling in the nonparametric bootstrap (Web Table S10)

Resampling subjects with replacement duplicates their event times, so every M4 bootstrap
resample contains tied event times. `timereg::aalen()` (the M4 arm of the main simulation)
breaks ties by adding a small random increment to every copy after the first and then
processes the copies one at a time. This experiment refits the **same B = 200 resamples of
the same canonical simulated datasets** with tied copies processed simultaneously, and
compares the two conventions resample by resample. It uses no true values and computes no
coverage.

The planned production run was stopped after 717 of 14,526 tasks (about 5% of its
replicates). Web Table S10 summarises the completed replicates (Panel A: production
replicates complete in every scenario of a sample size; Panel B: the pilot).

## Files

| File | Role |
|---|---|
| `tie_exact.R` | the two tie-exact refits: (E1) a from-scratch weighted least-squares Aalen fit with simultaneous ties; (E2) `timereg` on the resample with each duplicated subject collapsed to one row weighted by its multiplicity |
| `prepare_fast.R` | a base-R `ste_prepare()` that returns identical columns, for the compute budget |
| `m4_core.R` | one canonical replicate, refitted both ways on the same resamples (`paired_rep()`); modes `full` and `advance` |
| `validate_tie_exact.R` → `validation.log` | checks before the run: identical preparation, E1 = timereg on tie-free data, E1 = E2 on resamples, jitter live |
| `make_tasks.R` → `tasks_pilot.rds`, `tasks_prod.rds` | the fixed task lists (pilot: canonical replicates 999–1000; production: replicates 1..K) |
| `run_m4.R` | resumable runner, one `.rds` per task under `reps/<stage>/` |
| `extract_canonical.R` → `canon_boot.rds`, `canon_nfail.rds` | the canonical (timereg) bootstrap rows the tie-exact arm is paired with |
| `canonical_rows.R` → `canonical_rows.rds` | full-precision published M4 rows (Table 3 SE ratios) |
| `se_subset.R` → `se_subset.rds`, `se_subset.csv`, `se_subset_table.tex`, `se_subset.log` | **Web Table S10** and the projected Table 3 SE ratios (0.908 for 0.919 at n = 300; 0.839 for 0.847 at n = 1000) |
| `make_boottie_tex.py` + `boottie_appendix_template.tex` → `boottie_appendix.tex` | fills the Web Appendix D block (S10 rows and the projection figures) from `se_subset_table.tex` and `se_subset.log` |

## Reproduce Web Table S10

Extract the companion research-inputs archive linked in the repository README at
the repository root. It supplies `reps/pilot/` and `reps/prod/` here (837 replicate
files, 12 MB) and the per-replicate record of the
main simulation, `data/exact/results/replicates_all.rds` (91 MB), from which
`canon_boot.rds` (32 MB) is rebuilt. Then, from the repository root:

```sh
Rscript dev/m4_ties/extract_canonical.R              # -> dev/m4_ties/canon_boot.rds, canon_nfail.rds
Rscript dev/m4_ties/se_subset.R > dev/m4_ties/se_subset.log
python  dev/m4_ties/make_boottie_tex.py              # -> boottie_appendix.tex
```

`M4_DATA` points `se_subset.R` at another directory holding `canon_boot.rds` and `reps/`, and
`M4_OUT` redirects its outputs.

## Regenerate the replicate files

```sh
Rscript dev/m4_ties/validate_tie_exact.R             # optional pre-run checks
Rscript dev/m4_ties/run_m4.R pilot                   # M4_WORKERS (default 6) PSOCK workers
Rscript dev/m4_ties/run_m4.R prod
```

A replicate at n = 300 takes about half a minute in `advance` mode; the full production
design is about 175 CPU-hours, which is why it was stopped. The original run sourced a
frozen copy of the package code. The current package defaults reproduce those computations
(`dev/heart_ties/verify_default_path.R` reproduces canonical replicates to 1e-12), so use
the package code shipped in this release. `M4_SNAPSHOT` can point to a separately frozen
code directory when running further experiments.

The tie-exact refits ran under R 4.5.1 with `timereg` 2.0.7; the `timereg` arm is the main
simulation run (R 4.4.0 on the UBC ARC Sockeye cluster).
