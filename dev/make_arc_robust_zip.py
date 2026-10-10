#!/usr/bin/env python3
"""Build a self-contained zip for the robustness re-run on ARC.

Deliberately NOT PowerShell Compress-Archive: that writes BACKSLASH path separators, which Linux
unzip mangles into files literally named "simulation\\slurm\\robust.slurm". Uses zipfile with
explicit forward-slash arcnames, and normalises every text member to LF so the .slurm and .sh run
on Linux regardless of how git checked them out on Windows.

    python dev/make_arc_robust_zip.py          -> dev/steCI_robust_arc.zip
"""
import glob
import os
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "dev", "steCI_robust_arc.zip")
TOP = "steCI_robust_arc"

# Every text member is LF-normalised. steCI/R is included because robustness_misspec.R and
# robustness_tvfx.R source it by relative path; robustness_2confounder.R is self-contained but
# ships alongside. No scenarios.R -- verified none of the three reads DGM_FIXED or T_HOR.
TEXT = [
    "dev/robust_arc.R",             # the runner (the only thing the cluster executes)
    "dev/robust_collect.R",         # scoring; runs at HOME, shipped for reference
    "dev/robustness_2confounder.R",  # sourced for its definitions, NOT run
    "dev/robustness_misspec.R",
    "dev/robustness_tvfx.R",
    "simulation/slurm/robust.slurm",
    "simulation/slurm/submit_robust.sh",
] + sorted(
    os.path.relpath(p, ROOT).replace("\\", "/")
    for p in glob.glob(os.path.join(ROOT, "steCI", "R", "*.R"))
)

# NOTHING truth-related is bundled, and that is the design, not an omission: the cluster job
# computes no truth and reads no truth. See the README section below.
BINARY = []

README = """# Robustness re-run on ARC

Self-contained. No git needed; nothing outside this folder is read.

## Why

The three appendix robustness studies (`2confounder`, `misspec`, `tvfx`) score raw
coverage against ONE frozen Monte-Carlo draw of their own truth -- the same
single-frozen-draw truth defect as the primary DGM. The frozen draws sit 0.48-0.79 of a single-draw SD
from the exact values, a fixed displacement shared by every replicate that does not
average out.

The truths are now exact for two of the three. The reason this needs a CLUSTER RUN rather
than a recompute is separate and duller: **all three scripts saved only 9 summary rows and
threw away the per-replicate estimates**, so there is nothing left to re-score.

## The one design decision worth knowing

**This job computes NO truth and reads NO truth.** It emits one row per replicate --
`(gamma_hat, se_model, se_row, se_cl)` -- and stops. Scoring happens at home, in
`dev/robust_collect.R`. Three things follow:

1. `tvfx`'s truth is still underived, and that does not block this run at all.
2. Any FUTURE truth revision is a seconds-long re-score, never another cluster job.
3. The exact-truth files are deliberately not in this bundle.

Seeding is preserved from the original scripts (`seed = 1000*gi + r`), verified to
reproduce the original replicates to 0.000e+00, so every change in a reported coverage is
attributable to the truth alone rather than to a different draw.

## Run this

    unzip steCI_robust_arc.zip && cd steCI_robust_arc

    bash simulation/slurm/submit_robust.sh smoke     # 24 tasks, minutes -- ALWAYS FIRST
    squeue -u $USER

Read the two numbers the smoke prints (measured RSS, seconds per chunk), then:

    rm -f simulation/output/robust/rob_*_R2_*.rds     # clear smoke output
    bash simulation/slurm/submit_robust.sh full       # 357 tasks

Bring the per-replicate chunks home and score them there:

    tar czf robust_out.tgz simulation/output/robust

    # on the workstation, from the repo root:
    Rscript dev/robust_collect.R

## Profiles

| Profile | Tasks | What |
|---|---|---|
| `smoke` | 24 | 2 replicates per cell, all three studies. Output tagged `R2`, so it CANNOT collide with production. |
| `full` | 357 | The deliverable. 2confounder 15, misspec 171, tvfx 171. |
| `misspec` | 171 | Only the study with the largest displacement (0.79 draw SD at alpha0=0.40). |
| `tvfx` | 171 | Only tvfx. |
| `2conf` | 15 | Only the two-confounder study. |

Checkpointed: an existing chunk is skipped, so a re-submit resumes in seconds. Killed or
pre-empted tasks cost only their own chunk.

## Sizing

Per-replicate cost MEASURED locally: 0.84 / 2.50 / 30.4 s at n = 300 / 1000 / 5000. ARC
cores measured ~4x slower per core during the truth run, so chunks (250 / 100 / 10
replicates by n) target ~15-25 min per task. `--time=1:00:00` is ~2.5x headroom, because
the 4x factor is extrapolated from a different workload.

Peak R heap MEASURED at 229 MB for the worst cell (misspec n=5000, 19,585 stacked rows).
`--mem=8G` is deliberate headroom; the smoke profile reports the real cgroup peak, so cut
it if that comes back small.

## If something looks wrong

The runner sources each study's ORIGINAL script up to its driver loop, so the DGM and fit
code is byte-identical to what produced the published tables -- it is not reimplemented
here. If a study's `grid <- expand.grid` line ever moves, the runner FAILS a `stopifnot`
rather than sourcing the driver loop. That is intentional.
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
    print("members:")
    with zipfile.ZipFile(OUT) as z:
        for i in z.infolist():
            print(f"  {i.filename}")
            if "\\" in i.filename:
                raise SystemExit("BACKSLASH in arcname -- Linux unzip would mangle this")


if __name__ == "__main__":
    main()
