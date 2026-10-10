#!/usr/bin/env python3
"""Build a self-contained zip for the max.clust sensitivity on ARC.

Deliberately NOT PowerShell Compress-Archive: that writes BACKSLASH path separators, which Linux
unzip mangles. Uses zipfile with explicit forward-slash arcnames, and normalises every text member
to LF so the .slurm and .sh run on Linux regardless of how git checked them out on Windows.

    python dev/make_arc_maxclust_zip.py          -> dev/steCI_maxclust_arc.zip
"""
import glob
import os
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "dev", "steCI_maxclust_arc.zip")
TOP = "steCI_maxclust_arc"

# steCI/R is the shipped estimator (dat_sim + ste_prepare + ci_analytic, the byte-identical cluster
# fit); simulation/scenarios.R gives make_scenarios(), DGM_FIXED, T_HOR. NOTHING truth-related is
# bundled: the cluster computes no truth -- scoring is dev/maxclust_collect.R at
# home, which reads dev/truth_cache_exact.rds there.
TEXT = [
    "dev/maxclust_arc.R",              # the runner (the only thing the cluster executes)
    "dev/maxclust_collect.R",          # scoring; runs at HOME, shipped for reference
    "simulation/scenarios.R",          # make_scenarios(), DGM_FIXED, T_HOR
    "simulation/slurm/maxclust.slurm",
    "simulation/slurm/submit_maxclust.sh",
] + sorted(
    os.path.relpath(p, ROOT).replace("\\", "/")
    for p in glob.glob(os.path.join(ROOT, "steCI", "R", "*.R"))
)

README = """# max.clust sensitivity on ARC

Self-contained. No git needed; nothing outside this folder is read.

## Why

The published M3 cluster-robust SE for the constant hazard difference (E1) is
`timereg::aalen(..., clusters = id)` at timereg's **default `max.clust = 1000`**. When the number
of subject clusters exceeds 1000, timereg quantile-bins them into 1000 groups, so at **n = 5000**
(5000 clusters > 1000) distinct subjects share bins; at n = 300 and n = 1000 the clusters are
exact (1000 is not < 300/1000), so pinning `max.clust = NULL` changes nothing there and they are
NOT run. The binned "meat" is UNBIASED for the exact subject-level meat -- subjects are
independent, so the omitted cross-cluster term is mean-zero -- so the shift from pinning is noise
about zero, not bias. The paper states this on a **20-replicate** check (mean +0.37%, sd 1.63%,
t = 1.01 on 19 df, p = 0.33). This job tightens that bound: it re-fits the cluster SE **with and
without the pin** at every n = 5000 cell over many replicates, so the mean shift and its induced
coverage change are pinned down.

**This does not change the published numbers.** The paper keeps the timereg default (so the code
reproduces the tables); this run only quantifies what pinning WOULD do.

## The one design decision worth knowing

**This job computes NO truth and reads NO truth**. It emits one row per
replicate -- `(gamma_hat, se_default, se_pinned)` -- and stops. The default fit does NOT pass
`max.clust`, so it is byte-identical to `ci_analytic.R`'s published cluster fit; the pinned fit
adds `max.clust = NULL`. `max.clust` changes ONLY the robust variance, never `gamma` (asserted in
the runner). Scoring -- the relative SE shift, its t-test, the per-cell shifts, and the induced
E1 coverage change against the exact truths -- happens at home in `dev/maxclust_collect.R`, which
reads `dev/truth_cache_exact.rds` there. So a truth revision is a seconds-long re-score, never a
cluster job.

The estimator is NOT reimplemented: the runner calls `dat_sim()` + `ste_prepare()` + the shipped
cluster `aalen` fit from `steCI/R`.

## Version note & self-check

timereg's default `max.clust` has been **1000** across recent versions; the binning at n = 5000
(5000 clusters > 1000) is confirmed both locally (timereg 2.0.7: per-replicate shifts of order 1%,
mixed sign, sd ~1.6% -- exactly the paper's 20-rep picture) and on the ARC run that produced the
published tables. `dev/maxclust_collect.R` prints a BINNING SANITY block -- the resolved default
`max.clust`, the `nclust` range, and the fraction of replicates with a bit-identical (exactly-zero)
contrast -- and warns loudly if the contrast did not fire, so a timereg whose default differs
cannot be mistaken for a clean confirmation. Read that block first.

## Cells and sizes

All **27 n = 5000 scenarios** (`make_scenarios()` filtered to n = 5000), the only sample size
timereg bins. FULL = 100 replicates/cell, 10 reps/chunk -> 270 array tasks. That is ~2700
comparisons, so the pooled mean-shift SE is ~0.03% -- an order of magnitude tighter than the
20-replicate check.

## Run this

    unzip steCI_maxclust_arc.zip && cd steCI_maxclust_arc

    bash simulation/slurm/submit_maxclust.sh smoke     # 4 tasks (2 reps/cell) -- ALWAYS FIRST
    squeue -u $USER

Read the two numbers the smoke prints (peak RSS, seconds/fit), then:

    rm -f simulation/output/maxclust/mc_*_R2_*.rds       # clear smoke output
    bash simulation/slurm/submit_maxclust.sh full        # 270 tasks

Bring the per-replicate chunks home and score them there:

    tar czf maxclust_out.tgz simulation/output/maxclust
    # on the workstation, from the repo root:
    Rscript dev/maxclust_collect.R      # -> dev/maxclust_results.rds

## Sizing

Two aalen fits per replicate at n = 5000 (~11-15 s/fit average measured locally; the exact
`max.clust=NULL` fit builds a larger meat and is the slower of the two, ~20 s). A 10-replicate task
is ~20 fits (~5 min), so `--time=0:30:00` is generous -- but the SMOKE measures it; read it first.
`--mem=4G` is ample (n = 5000 is well under the n = 1e4 plim run's ~630 MB peak). Checkpointed on
file existence: a re-submit skips completed chunks in seconds; a killed task costs only its own
chunk (atomic temp+rename write).
"""


def lf(path):
    with open(path, "rb") as fh:
        return fh.read().replace(b"\r\n", b"\n")


def main():
    missing = [m for m in TEXT if not os.path.exists(os.path.join(ROOT, m))]
    if missing:
        raise SystemExit("missing members:\n  " + "\n  ".join(missing))

    with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as z:
        for m in TEXT:
            z.writestr(f"{TOP}/{m}", lf(os.path.join(ROOT, m)))
        z.writestr(f"{TOP}/README.md", README.encode("utf-8"))

    print(f"wrote {OUT}  ({len(TEXT) + 1} members, {os.path.getsize(OUT) / 1024:.0f} KB)")
    with zipfile.ZipFile(OUT) as z:
        for i in z.infolist():
            print(f"  {i.filename}")
            if "\\" in i.filename:
                raise SystemExit("BACKSLASH in arcname -- Linux unzip would mangle this")


if __name__ == "__main__":
    main()
