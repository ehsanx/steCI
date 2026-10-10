#!/usr/bin/env python3
"""Build a self-contained zip for the pooled-plim empirical run on ARC.

Deliberately NOT PowerShell Compress-Archive: that writes BACKSLASH path separators, which Linux
unzip mangles into files literally named "simulation\\slurm\\plim.slurm". Uses zipfile with
explicit forward-slash arcnames, and normalises every text member to LF so the .slurm and .sh run
on Linux regardless of how git checked them out on Windows.

    python dev/make_arc_plim_zip.py          -> dev/steCI_plim_arc.zip
"""
import glob
import os
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "dev", "steCI_plim_arc.zip")
TOP = "steCI_plim_arc"

# Every text member is LF-normalised. steCI/R is the shipped estimator (dat_sim + ste_estimate);
# simulation/scenarios.R gives DGM_FIXED and T_HOR (the two repo copies are byte-identical, so the
# shallow one ships and dev/plim_arc.R resolves it). NOTHING truth-related is bundled: the cluster
# computes no truth -- scoring is dev/plim_collect.R, run at home.
TEXT = [
    "dev/plim_arc.R",              # the runner (the only thing the cluster executes)
    "dev/plim_collect.R",          # scoring; runs at HOME, shipped for reference
    "simulation/scenarios.R",      # DGM_FIXED, T_HOR
    "simulation/slurm/plim.slurm",
    "simulation/slurm/submit_plim.sh",
] + sorted(
    os.path.relpath(p, ROOT).replace("\\", "/")
    for p in glob.glob(os.path.join(ROOT, "steCI", "R", "*.R"))
)

BINARY = []

README = """# Pooled-plim empirical run on ARC

Self-contained. No git needed; nothing outside this folder is read.

## Why

This run checks whether the residual long-horizon E2 point-estimate bias is a
NON-VANISHING pooled-origin offset (the additive MSM pools the cumulative coefficients B(.)
across trial origins, which the DGM mildly violates) or a finite-sample effect that shrinks
with n. The n -> infinity value is now settled by `dev/plim_closed_form.R` on the workstation (it
solves the population normal equations from accumulated moments, Rao-Blackwellised, and never runs
an aalen fit): the asymptotic bias is small (~0.007 at tau=5 in the worst cell) and is NOT
origin-pooling. This ARC job does ONE complementary thing: an INDEPENDENT-CODE check that the
SHIPPED `ste_estimate` (the real `timereg::aalen` path, no Rao-Blackwellisation) gives a clearly
POSITIVE finite-n bias at n=1e4 (~0.024 +/- 0.006), consistent with that positive asymptote --
closing the "is your closed form the same object as your package?" gap.

**n is capped at 1e4.** `timereg::aalen` is empirically ~O(n^2) here (85-122 s/fit at n=1e4, 985
s/fit at n=3e4, ~3 h at n=1e5), and -- decisively -- the n>=3e4 arm carries NO resolvable signal:
the finite-n bias moves only ~0.002 between n=5000 and n=3e4 while the tau=5 SE is ~0.008, so it
cannot discriminate the hypotheses. Do not re-add higher n.

## The one design decision worth knowing

**This job computes NO truth and reads NO truth.** It emits one row per replicate --
`(gamma_hat, mrd_hat[1..5])` from the shipped `ste_estimate()` -- and stops. Scoring against the
exact closed-form truths, the vanish-vs-plateau read, and the closed-form cross-check all happen
at home in `dev/plim_collect.R`. So any future truth revision is a seconds-long re-score, never
another cluster job, and the exact-truth files are deliberately not in this bundle.

The estimator is NOT reimplemented: the runner calls `dat_sim()` + `ste_estimate()` from the
shipped `steCI/R`, byte-identical to the code behind the published tables. A disagreement with
the closed form is therefore a derivation issue, not a transcription slip.

## Cells and sizes

Four cells (`dev/plim_arc.R`), all at prevalence-high (gamma.0 = 0):
  1 evlow_cfstrong_pvhigh  (alpha.0 0.15, gamma.L 3.0) -- the worst observed-bias cell
  2 evlow_cfweak_pvhigh    (alpha.0 0.15, gamma.L 0.5) -- ~zero-bias control
  3 evlow_cfmod_pvhigh     (alpha.0 0.15, gamma.L 1.5) -- fills the confounding axis
  4 evhigh_cfstrong_pvhigh (alpha.0 0.40, gamma.L 3.0) -- event-rate check
n = 10000 only, 200 replicates/cell (SE ~0.006 on the tau=5 estimate). The E2
truth depends on alpha.0 ONLY.

## Run this

    unzip steCI_plim_arc.zip && cd steCI_plim_arc

    bash simulation/slurm/submit_plim.sh smoke     # 4 tasks (2 reps/cell) -- ALWAYS FIRST
    squeue -u $USER

Read the two numbers the smoke prints (peak RSS, seconds/fit), then:

    rm -f simulation/output/plim/plim_*_R2_*.rds       # clear smoke output
    bash simulation/slurm/submit_plim.sh full          # 160 tasks

Bring the per-replicate chunks home and score them there:

    tar czf plim_out.tgz simulation/output/plim
    # on the workstation, from the repo root:
    Rscript dev/plim_collect.R

## Sizing -- MEASURED on ARC (job 12342697, 2026-07-22), not extrapolated

Per-fit wall cost **85-122 s at n=1e4**, peak RSS **470-631 MB**, hence `--mem=4G`. 160 tasks x 5
replicates/chunk = ~10 min/task, 0:30:00 walltime, ~22 CPU-hr total.

The first version used 15-25 replicates/task and n up to 1e5, sized off an `n^0.6` extrapolation
from a local 5000->10000 probe; the true scaling is ~O(n^2), so those tasks hit the walltime, and
the n>=3e4 arm carried no resolvable signal anyway. Now: n=1e4 only, many short tasks (they
backfill). Do NOT raise n or `chunk_reps()`.

Checkpointed on file existence: a re-submit skips completed chunks in seconds; a killed task
costs only its own chunk (atomic temp+rename write).
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
        for src, dst in BINARY:
            with open(os.path.join(ROOT, src), "rb") as fh:
                z.writestr(f"{TOP}/{dst}", fh.read())
        z.writestr(f"{TOP}/README.md", README.encode("utf-8"))

    n = len(TEXT) + len(BINARY) + 1
    print(f"wrote {OUT}  ({n} members, {os.path.getsize(OUT) / 1024:.0f} KB)")
    with zipfile.ZipFile(OUT) as z:
        for i in z.infolist():
            print(f"  {i.filename}")
            if "\\" in i.filename:
                raise SystemExit("BACKSLASH in arcname -- Linux unzip would mangle this")


if __name__ == "__main__":
    main()
