# Handoff: benchmark session (2026-08-11/12, branch `speed-up3`, PR #43)

Started as "benchmark the RNG speed-up against the field and put it in the PR".
Ended up finding three bugs, two of them not in this branch. Full results:
`2026-08-11-ransac-benchmark.md`. This document is the state of play and what
is still open.

## What landed

| commit | what |
|---|---|
| `d4e7ef4` | design spec (`docs/superpowers/specs/2026-08-11-ransac-benchmark-design.md`) |
| `a8caab8` | `benchmarks/` harness: self-contained F + H accuracy-vs-compute benchmark |
| `feb9e24` | first full results + report |
| `0ae6b99` | figures-first write-ups; `rng_cost.c` |
| `e8c028a` | corrected the macOS compiler guidance (prefer Clang, not GCC) |
| `899af7f` | PROSAC score-orientation fix; published-release comparison |
| `f349a6c` | **macOS LAPACK fix**; CI to manylinux_2_28 + OpenBLAS; `PYDEGENSAC_CMAKE_ARGS` |
| `6d1b1ba` | perf profile of the wheel-vs-local gap |

## Headline results (Linux, x86-64, glibc)

- Branch vs master: **1.19x (F), 1.20x (H)** aggregate, 1.22x / 1.39x at the
  largest budgets. No accuracy difference is significant at any budget.
- The **4.2x / 2.2x** in the PR title are macOS golden-pair figures measured on
  a build where LAPACK was compiled out entirely (below). They describe a
  configuration no Linux or Windows user has run.
- Against the field, pydegensac is measurably behind poselib-prosac on both
  problems: -0.023 mAA on F at 1.5x the cost, -0.012 on H at 5x
  `cv2.USAC_MAGSAC`'s cost.

## Three bugs found

1. **macOS never called LAPACK** (fixed in `f349a6c`). `lapwrap.c` guarded both
   the prototypes and the calls behind `#ifdef _WIN32 || __linux__`; macOS
   defines neither. Every least-squares F refit returned the identity matrix
   (`singulF`'s error path), and every H least-squares refit copied an
   untouched covariance matrix. LO-RANSAC's local optimisation has been dead on
   macOS. Costs 0.11 mAA on F and 0.21 on H.
2. **`pydegensac==0.1.2` silently returns every point as an inlier under
   numpy 2.x** — the old-pybind11 bug that `0.2` was yanked for. 0.1.2 is *not*
   yanked and is what anything installed before 2026 has pinned. Not fixable in
   this repo; needs a yank on PyPI.
3. **PROSAC was sampling worst-first on HPatchesSeq** (fixed in `899af7f`) —
   benchmark-only bug. The tutorial archives disagree about which way
   `match_conf` points; `check_scores.py` measures it against ground truth.

## Open decisions (maintainer's call, deliberately not made)

- **PR title** still claims 4.2x / 2.2x. Those numbers are real only for the
  macOS no-LAPACK build.
- **macOS golden baselines** encode the broken behaviour. Every golden test on
  macOS will fail until regenerated, and the regeneration should be reviewed by
  a human rather than rubber-stamped — the new outputs should be *better*, and
  that is worth confirming pair by pair.
- **glibc floor**: moving Linux wheels to manylinux_2_28 drops CentOS 7 /
  Ubuntu 18.04. Reversible by deleting one line
  (`CIBW_MANYLINUX_X86_64_IMAGE`).
- **cp38**: the main Linux job still uses cibuildwheel 2.21.2 because 3.x
  dropped cp38. If cp38 support is dropped, the two Linux jobs can merge and
  the `numpy<2.3` workarounds go away.

## Open threads

### 1. The unexplained 75% of the wheel gap (next experiment, needs Docker)

Published Linux wheels are **1.10x (F) / 1.30x (H) slower than the same source
built locally**. Reproduced interleaved, three times each. perf says the extra
time is inside pydegensac's own compiled code (0.456 s vs 0.312 s isolated,
1.46x), not libc and not LAPACK.

Ruled out by direct experiment — **do not re-test these**: LAPACK
implementation (2%), optimisation level (`-O2` 0.44 s vs `-O3` 0.43 s), GCC
major version (conda gcc 10.4 ≈ gcc 13.3), symbol visibility (wheel exports 425
symbols vs 121 locally, but hiding them buys 5%), assertions, hardening flags,
instruction mix (59,915 vs 60,779, both all-scalar).

Not yet isolated: the wheel records `GCC 10.2.1`, RedHat's **devtoolset** build,
while the gcc-10 control was conda's 10.4.0 — a different GCC with different
defaults and a CentOS 7 baseline. The experiment, once Docker is available:

```bash
# in each of quay.io/pypa/manylinux2014_x86_64 and manylinux_2_28_x86_64:
#   build tag v_0.2.2, copy the .so out, time it with the isolated driver
#   against the local gcc 13.3 build, threads pinned
```

If manylinux2014 reproduces 0.71 s and manylinux_2_28 does not, the image move
in `f349a6c` is worth ~30% to every Linux `pip install` user rather than the
~8% currently claimed, and the report and PR should be updated to say so. If
both are slow, the image move buys only the LAPACK share and something else is
going on.

### 2. macOS/M1 re-benchmark

An M1 run was started against `899af7f`, i.e. *before* the LAPACK fix — those
numbers describe the crippled build and should be discarded or relabelled. Worth
re-running after `f349a6c`. Check which build you have with:

```bash
nm -u $(python -c "import pydegensac,glob,os;print(glob.glob(os.path.dirname(pydegensac.__file__)+'/*.so')[0])") | grep -E 'dgesvd|dsyev'
```

Empty output = the crippled build. Also worth running `rng_cost` there: it
should report a per-iteration saving several times the 406 ns measured on
glibc, which is the claim that the 4.2x/2.2x is a macOS-libc artefact.

### 3. Smaller items

- `lapwrap.c`'s workspace memoisation measured **neutral** (18.4 ms before and
  after on H). Kept for tidiness. If someone later profiles LAPACK as hot on
  another platform, this is already done.
- The F ratio grid tops out at 1.0 and the H one at 64 px; both optima are now
  interior, but `cv2-magsac` on HPatchesSeq was still creeping upward at 64 px
  (0.9324 -> 0.9331). Immaterial, but it is the one boundary left.

## Working with the harness

`benchmarks/README.md` is the reference. Things a fresh session will trip over:

- **`benchmarks/.ab/` holds ~25 throwaway build trees** (one per version /
  compiler / flag variant) plus three conda prefixes and three git worktrees.
  All gitignored. `rm -rf benchmarks/.ab` to reclaim; `run_ab.sh` and
  `run_releases.sh` rebuild what they need. `benchmarks/data` is another 1.5 GB.
  Two git worktrees are still registered (`worktree-base`, `wt-local-0.2.2`,
  `wt-local-master`) — `git worktree remove --force` them if `.ab` is deleted.
- **This machine cannot build pydegensac from a plain `pip install .`** without
  the conda prefix at `benchmarks/.ab/env` — there was no CMake and no LAPACK
  until 2026-08-12, when `pkg-config`, `cmake`, `liblapack-dev` and
  `libopenblas-dev` were installed system-wide. A plain build now works.
- **Golden tests fail in exact mode on Linux** (20 of 33) — expected, the
  baselines are macOS-captured. Use `PYDEGENSAC_GOLDEN_EXACT=0`, which is what
  CI does. To check a C change is output-preserving on Linux, compare *seeded*
  outputs instead; that is how `f349a6c` was verified byte-identical.
- **Two measurement traps**, both of which produced wrong conclusions before
  being caught:
  - profiling `run.py` shows numpy + cv2 at 85% — that is the reprojection
    metric, not the estimator. Profile an estimator-only driver.
  - unpinned OpenBLAS burns 20-43% in `__sched_yield` spinning on 9x9
    matrices. `run.py` pins threads; ad-hoc scripts must too.
- **perf on this box**: built from the WSL2 kernel tree, needs
  `LD_LIBRARY_PATH=$HOME/miniconda3/envs/py313/lib` (it linked conda's
  libpython). WSL2 exposes **no hardware PMU** — software events only, so no
  instruction counts or IPC.
- **Statistics**: mAA differences below ~0.02 (F) / ~0.037 (HPatchesSeq) are
  run-to-run scatter. Use the paired bootstrap in `report.py`; the marginal CIs
  are far wider and will tell you nothing separates.
- **EVD (8 test pairs) cannot rank anything.** It is reported because ds-sac's
  protocol includes it.
