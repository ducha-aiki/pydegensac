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

- ~~**PR title** still claims 4.2x / 2.2x~~ — retitled 2026-08-12 to lead with
  the LAPACK fix and the measured 1.2-1.4x.
- ~~**macOS golden baselines** encode the broken behaviour~~ — regenerated on
  M1 against the fixed build (`edabc68`); the gate is green in exact mode
  (33 passed). Reviewed pair by pair as asked: inlier counts rise on 17 of 18
  (seed, pair) cases — 99 -> 128 on `adam`, 61 -> 115 on `face`, 60 -> 71 on
  the first reichstag pair — with GT error flat or better on H and within
  0.02 px on F. The one exception (`05534141-05545431`, seed 2024) is
  357 -> 354 inliers and +0.001 px, i.e. noise. `05466646-05534141` no longer
  passes capture sanity with LO restored and is replaced by
  `06373813-06639257`.
- **glibc floor**: moving Linux wheels to manylinux_2_28 drops CentOS 7 /
  Ubuntu 18.04. Reversible by deleting one line
  (`CIBW_MANYLINUX_X86_64_IMAGE`).
- **cp38**: the main Linux job still uses cibuildwheel 2.21.2 because 3.x
  dropped cp38. If cp38 support is dropped, the two Linux jobs can merge and
  the `numpy<2.3` workarounds go away.

## Open threads

### 1. The wheel gap — CLOSED 2026-08-12, it is `pinvJ`

Docker was available after all. Ran the experiment: rebuilt the same source in
both manylinux images, copied the bare un-repaired `.so` out, timed all arms
interleaved on one P-core against the host's OpenBLAS so LAPACK is constant.

**manylinux2014 reproduces the wheel; manylinux_2_28 does not.** The image move
in `f349a6c` is worth **~1.25x on H and ~1.10x on F** to every Linux `pip
install` user — essentially the whole gap, not the ~8% previously claimed.
Nothing is left unexplained: manylinux2014 + the wheel's own vendored LAPACK
reproduces the published wheel to within 1%.

Root cause is a single function. GCC 10.2.1 (CentOS 7 devtoolset) spills
`pinvJ`'s `pJ[i] /= N` loop to memory and serialises its four `divpd` behind
load/store round-trips; GCC 14 keeps it in registers and overlaps them.
`pinvJ` is 2.7x slower and accounts for the *entire* process-level difference
(926 of 921 net samples). It lives in `Htools.c` and is called only from there,
which independently explains why H loses 30% and F only 5%.

Also rejected, on top of the earlier list: symbol visibility, despite the 425
vs 125 exported-symbol correlation (`-fvisibility=hidden
-fno-semantic-interposition` in manylinux2014 buys 1.6%).

Full write-up, tables and disassembly: "The wheel gap is `pinvJ`" in
`2026-08-11-ransac-benchmark.md`. The harness is checked in at
`benchmarks/toolchain/` with its own README — `build_in_image.sh` builds a ref
inside an image, `iso_bench.py` is the estimator-only fixed-seed driver,
`run_arms.sh` interleaves arms, `summarize.py` prints medians and verifies the
arms agree on inliers. Build outputs land in the disposable `benchmarks/.ab/`.

**Fell out of this, and it matters for CI:** `v_0.2.2` and `master` *cannot
build* in manylinux_2_28 — `Ftools.c` calls `dgeqp3_` with no prototype and
gcc 14 makes that an error. The `lapwrap.h` declaration on this branch is what
makes the image move viable, so `f349a6c`'s two halves have to ship together.
The macOS half of the same gap (declaration and calls both `#ifdef`-guarded
away, so `dgeqp3_` never ran there either) was closed independently in
`614e2e0` while this was being measured.

### 2. macOS/M1 re-benchmark — DONE (`6b11a8f`)

Re-run on M1 against `f349a6c`, three arms (`master`, `master` + the fixed
`lapwrap.c`, branch), full roster on both problems. Results in the M1 section
of `2026-08-11-ransac-benchmark.md`. Summary:

- **1.35x (F) / 1.18x (H HPatchesSeq) at equal correctness** — that is the
  macOS number, replacing 4.2x / 2.2x.
- The pre-fix M1 run (against `899af7f`) reproduced the crippled build natively
  and matched the Linux `-U__linux__` simulation closely: F best mAA 0.3133 vs
  0.3272 simulated, H 0.6924 vs 0.7021. The bug cost macOS 0.11 / 0.21 mAA.
- `rng_cost` on M1: **1093 ns saved per F iteration, 1079 ns per H**, against
  406 / 366 on glibc. 2.7x, so the libc lock is confirmed as the platform
  difference and ARM is not a factor.

Check which build you have with:

```bash
nm -u $(python -c "import pydegensac,glob,os;print(glob.glob(os.path.dirname(pydegensac.__file__)+'/*.so')[0])") | grep -E 'dgesvd|dsyev'
```

Empty output = the crippled build.

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
