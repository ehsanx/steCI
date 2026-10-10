#!/usr/bin/env python3
"""Build a self-contained zip for the M8 (E1 target) ARC run.

Deliberately NOT PowerShell Compress-Archive: that writes BACKSLASH path separators, which Linux
unzip mangles into files literally named "simulation\\slurm\\truth.slurm". Uses zipfile with explicit
forward-slash arcnames, and normalises every text member to LF so the .slurm and .sh run on Linux
regardless of how git checked them out on Windows.

    python dev/make_arc_truth_zip.py          -> dev/steCI_truth_arc.zip
"""
import io
import os
import sys
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "dev", "steCI_truth_arc.zip")
TOP = "steCI_truth_arc"

# Text members are LF-normalised; binary members are copied byte-for-byte.
TEXT = [
    "dev/truth_arc.R",
    "dev/truth_collect.R",
    "simulation/slurm/truth.slurm",
    "simulation/slurm/submit_truth.sh",
    "simulation/scenarios.R",
]
BINARY = [
    # 749 bytes. Lets truth_collect.R print the "vs cached" column on ARC.
    ("data/truth_cache.rds",
     "dev/truth_cache_canonical.rds"),
]

README = """# M8 -- re-estimate the E1 target `gamma_true` on ARC

Self-contained. No git needed; nothing outside this folder is read.

## Why

The shipped E1 target is ONE frozen Monte-Carlo draw. Its replication SD is
0.00254 / 0.00382 / 0.00629 at alpha.0 = 0.15 / 0.25 / 0.4, against an E1 bias of
about 0.0057 -- so the noise in the target is comparable to the quantity it
measures, and being one draw shared across every replicate and scenario it is a
fixed displacement that does not average out.

## Run this

    unzip steCI_truth_arc.zip && cd steCI_truth_arc
    bash simulation/slurm/submit_truth.sh quick      # 189 tasks, ~12-20 min each
    squeue -u $USER

    module load gcc/9.4.0 r/4.4.0
    export R_LIBS=$HOME/R/x86_64-pc-linux-gnu-library/4.4
    Rscript dev/truth_collect.R

`quick` averages K=2000 draws at n_hd=3e4 and is the whole point of the exercise:
it takes the target's Monte-Carlo SE from ~110% of the E1 bias down to ~2.5%.
Bring back `simulation/output/truth/`, or just the printed table.

## Profiles

| profile | tasks | K    | n_hd            | time    | what it is for |
|---------|-------|------|-----------------|---------|----------------|
| quick   | 189   | 2000 | 3e4             | 1:00:00 | THE deliverable: the variance fix |
| scale   | 9     | 2    | 3e4, 4.5e4, 6e4 | 2:00:00 | measure cost/memory growth in n_hd |
| smoke   | 3     | 2    | 3e4             | 1:00:00 | plumbing check after editing scripts |
| full    | 567   | 2000 | 3e4, 1e5, 3e5   | 3:00:00 | the sweep -- DO NOT run unsized |

Run `full` only after reading `scale`. At ARC's measured speed the 3e5 arm is
plausibly ~19 h per task and well beyond --mem=64G.

## Measured on ARC, 2026-07-19

Per-draw wall cost at n_hd=3e4: 372 / 552 / 697 s at alpha.0 = 0.15 / 0.25 / 0.4.
The same draws cost 88 / 130 / 165 s on a laptop, so ARC cores are ~4x slower per
core -- walltimes here are sized off the ARC numbers. Peak RSS 628-985 MB per
worker, i.e. ~20-32 GB at 32 workers, so 64G is comfortable at n_hd=3e4. Memory
above 3e4 is unmeasured; that is what `scale` is for.

## Gotchas that already bit us

- `TRUTH_NHD` is **colon**-separated. `sbatch --export` is comma-separated, so a
  comma inside a value silently truncates the list -- on the first smoke run the
  sweep never ran and 6 of 9 tasks were no-ops, with no error.
- `TRUTH_CHUNK` must stay <= `--cpus-per-task`. CHUNK > cores is the bug that
  walltime-killed 10 of 81 cells in the main simulation.
- Output filenames embed alpha.0, n_hd, K and CHUNK, so a smoke chunk can never be
  silently reused as production. Everything is checkpointed: re-running a profile
  skips completed chunks in seconds.
- Needs R packages `timereg` and `survival`.
"""


def main():
    missing = [p for p in TEXT if not os.path.exists(os.path.join(ROOT, p))]
    missing += [s for s, _ in BINARY if not os.path.exists(os.path.join(ROOT, s))]
    if missing:
        sys.exit("missing inputs:\n  " + "\n  ".join(missing))

    steci = sorted(
        "steCI/R/" + f
        for f in os.listdir(os.path.join(ROOT, "steCI", "R"))
        if f.endswith(".R")
    )
    if not steci:
        sys.exit("no steCI/R/*.R found")

    n_text = 0
    with zipfile.ZipFile(OUT, "w", zipfile.ZIP_DEFLATED) as z:
        for rel in TEXT + steci:
            raw = io.open(os.path.join(ROOT, rel), "rb").read()
            data = raw.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
            z.writestr(TOP + "/" + rel, data)          # forward slashes, always
            n_text += 1
        for src, dst in BINARY:
            z.write(os.path.join(ROOT, src), TOP + "/" + dst)
        z.writestr(TOP + "/README_ARC.md", README.encode("utf-8"))

    # Verify what we actually produced, rather than trusting the writer.
    with zipfile.ZipFile(OUT) as z:
        names = z.namelist()
        bad_sep = [n for n in names if "\\" in n]
        crlf = []
        for n in names:
            if n.endswith((".R", ".sh", ".slurm", ".md")):
                if b"\r" in z.read(n):
                    crlf.append(n)
        print("wrote %s  (%.1f KB, %d members)" % (OUT, os.path.getsize(OUT) / 1024.0, len(names)))
        print("  text members LF-normalised: %d" % n_text)
        print("  backslash separators: %d %s" % (len(bad_sep), "OK" if not bad_sep else bad_sep))
        print("  members still containing CR: %d %s" % (len(crlf), "OK" if not crlf else crlf))
        need = {TOP + "/dev/truth_arc.R", TOP + "/simulation/slurm/truth.slurm",
                TOP + "/simulation/slurm/submit_truth.sh", TOP + "/simulation/scenarios.R",
                TOP + "/dev/truth_collect.R", TOP + "/README_ARC.md"}
        miss = need - set(names)
        print("  required members present: %s" % ("OK" if not miss else miss))
        print("  steCI sources: %d" % len(steci))
        if bad_sep or crlf or miss:
            sys.exit("ZIP FAILED VERIFICATION")


if __name__ == "__main__":
    main()
