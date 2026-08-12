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
- ~~**glibc floor**: moving Linux wheels to manylinux_2_28 drops CentOS 7 /
  Ubuntu 18.04~~ — **decided 2026-08-12: keep the floor.** The main Linux job
  is back on the default manylinux2014, so CentOS 7 / Ubuntu 18.04 users can
  still install. Sections 2 and 2b measured the image choice as
  performance-neutral on both counts — the compiler term died with `pinvJ`'s
  array, and the LAPACK term died with the `openblas-devel` line rather than
  the image — so the floor now costs nothing. The cp314 job stays on
  manylinux_2_28 because numpy publishes no manylinux2014 wheel for cp314.
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

## For the Linux side: what landed on M1 and what needs checking there

A second optimisation pass landed on 2026-08-12 after the LAPACK fix, because
the profile that guided the original RNG work had been taken on the build where
LO was dead and pointed at the wrong functions. On M1 the series is worth
**F 1.88x, H 2.49x** (19.2 -> 36 and 250 -> 621 estimator calls per second).
Commits `8648752` through `15445ef`.

Everything below is measured on M1 only. Four things want Linux eyes, in
priority order.

### 1. LTO under the wheel toolchains (`d9cb2a4`) -- ANSWERED: it was broken (`01d899b`)

**Every LTO build on Linux failed to import**, with `undefined symbol: mattr` --
local gcc 13.3 and both manylinux images; every non-LTO build was fine. Not the
suspected `pinvJ` miscompile, and not image-specific: a link failure. Since LTO
is on by default wherever `check_ipo_supported()` passes, and it passes in both
images, the next wheel matrix would have published broken wheels for every
Linux Python.

Latent bug that LTO exposed, not a compiler problem. `mattr` is defined in
`matutls` and called from `pydegensac_support`, but `target_link_libraries`
listed `matutls` first, and a static archive only yields the members needed by
what the linker has already seen. Without LTO the flat order happened to work;
GCC's linker plugin resolves archive members differently and does not. Fixed by
declaring the dependency (`pydegensac_support PUBLIC matutls`) so CMake orders
and repeats the archives itself.

With LTO on after the fix: imports in all three toolchains, 33 tests pass, H
inlier checksum identical to the non-LTO build (142796), ~3% on H locally. So
LTO stays on and is bit-exact here too. `PYDEGENSAC_NO_LTO=1` is what bisected
it -- worth keeping.

Original note follows.



CMakeLists now enables `CMAKE_INTERPROCEDURAL_OPTIMIZATION` when
`check_ipo_supported()` says yes, with `-fno-strict-aliasing` alongside it.
Worth 5% on F here and bit-exact.

**The risk is specific**: LTO under manylinux2014's devtoolset gcc 10.2.1 is
exactly the configuration `a3b4d4e` found miscompiling `pinvJ`. Build the wheel
matrix, and if anything looks wrong set `PYDEGENSAC_NO_LTO=1` to bisect --
that switch exists for this. The `-fno-strict-aliasing` is load-bearing, not
decoration: this C core predates the rule and the C++ half already suppresses
the warning.

### 2. Does the wheel gap survive `pinvJ_sqsum`? (`8648752`)

`a3b4d4e` pinned the entire published-wheel slowdown on gcc 10.2.1 spilling
`pinvJ`'s `pJ[i] /= N` loop to memory and serialising four `divpd`. **That loop
and that array are gone** -- the eight divisions factor into one, and nothing
is written to memory to be spilled.

Re-run `benchmarks/toolchain/run_arms.sh` against `15445ef`. If the two images
converge, the manylinux_2_28 move becomes optional rather than load-bearing,
and the CentOS 7 / Ubuntu 18.04 glibc floor can be kept. That is a real
decision that this change may have taken off the table.

**ANSWERED: they converge.** Same protocol as before, bare `.so` against one
host OpenBLAS, LTO on, at `01d899b`:

| arm | H ms/pair | | F ms/pair | |
|---|---|---|---|---|
| local, gcc 13.3 | 3.24 | 1.00x | 15.34 | 1.00x |
| manylinux2014, gcc 10.2.1 | 3.55 | 1.095x | 15.24 | 0.99x |
| manylinux_2_28, gcc 14.2.1 | 3.54 | 1.092x | 15.36 | 1.00x |

The two images are now indistinguishable from each other (3.550 vs 3.542 on H)
where the old code had them at 1.41x vs 1.055x. Removing the `pJ[]` array
removed what gcc 10.2.1 was spilling, so the `pinvJ` finding is history rather
than a live constraint.

**This measures the compiler, not the wheel.** It deliberately excludes what
`auditwheel` vendors, which is the other half of the image choice -- see the
vendored-BLAS section below before trading `manylinux_2_28` away for the glibc
floor.

Unexplained leftover, small and image-independent: both containers sit ~9%
above local gcc 13.3 on H while matching each other exactly, despite being
gcc 10.2 and gcc 14.2. Two very different compilers agreeing with each other
and differing from a third points at something environmental rather than
codegen. Not chased.

### 2b. The other half: what `auditwheel` vendors -- also converged

Built real repaired wheels in both images with
`benchmarks/toolchain/build_wheel_in_image.sh` (CI's `CIBW_BEFORE_ALL_LINUX`
verbatim, then `auditwheel repair`), installed them and timed them as shipped,
so each uses its own vendored BLAS rather than the host's:

| arm | H ms/pair | | F ms/pair | |
|---|---|---|---|---|
| manylinux2014 wheel, vendors OpenBLAS 0.3.3 | 3.642 | 1.000 | 16.54 | 1.000 |
| manylinux_2_28 wheel, vendors OpenBLAS 0.3.15 | 3.623 | 0.995 | 16.55 | 1.001 |
| local bare `.so` + host OpenBLAS | 3.260 | 0.895 | 15.20 | 0.919 |

**Indistinguishable** -- 0.5% on H and 0.1% on F, both inside the run-to-run
spread. Three years of OpenBLAS buys nothing here, which is what you would
expect when every LAPACK call is on a 9x9 matrix, far below the size where
kernel work pays.

Worth noting *why* this is now a comparison of two OpenBLAS versions rather
than OpenBLAS against the 2012 reference build: `f349a6c` added
`openblas-devel` to `CIBW_BEFORE_ALL_LINUX`, so manylinux2014 no longer vendors
reference LAPACK 3.4.2. That single line, not the image move, is what killed
the LAPACK term.

**So the image choice is now performance-neutral on both counts** -- compiler
(2 above) and vendored libraries (here). Whether to keep `manylinux_2_28` is a
pure compatibility and tooling question now: going back to manylinux2014 would
restore the CentOS 7 / Ubuntu 18.04 floor at no measurable speed cost, but it
still drags the cibuildwheel-2.x/cp38 split and the numpy `manylinux_2_17`
situation along with it, which is what the workflow comments already weigh.
Nothing here forces the decision either way; it just removes speed as an
argument.

The ~9-10% wheel-vs-local residual shows up identically in both wheels, so it
is not the vendored BLAS either. Same unexplained environmental term as above.

### 3. Re-measure the series on Linux -- DONE, and the M1/Linux gap is explained

`08464ca` vs `6e72f7f`, full roster, then the pydegensac arms again to separate
signal from scatter. **F 2.18x / 2.20x, H HPatchesSeq 1.73x / 1.72x** aggregate;
timing reproduces to within 1% at every budget. Report tables and figures
updated in place (`results/time_maa_{f,h}.png` were regenerated by the run, so
the old 1.19x/1.20x text next to them had to go).

They were right that the numbers do not transfer, but the reason is mostly the
denominator, not the optimisations: on F the branch lands at the *same absolute
cost* on both platforms (144.6 ms Linux, 144.8 M1) and all of the ratio
difference is M1's slower `master`. On H, M1's `master` is 1.44x slower than
Linux's while M1's branch is 1.17x faster, and 1.44 x 1.17 = 1.68 = the
2.91/1.73 gap -- so ~2/3 of it is macOS having started further behind, via the
`srandom()` lock that H feels most (1079 ns/iteration saved on M1 vs 366 on
glibc, and H runs many cheap iterations). Written up under "Why M1 reads 2.91x
on H where Linux reads 1.73x".

**Accuracy: neutral, but only the repeat proves it.** Run 1 threw two `*`
marks, both favouring base -- F@1000 (-0.0209) and H@400 (-0.0404) -- and F's
low budgets were five consecutive negatives, which is a ~6% coincidence and
worth checking rather than dismissing. The repeat disagrees completely: F@1000
comes back -0.0016, H@400 comes back **+0.0264**, nothing is significant
anywhere on either problem, and the sign pattern is different. M1 shows no
low-budget deficit either. If you quote "no accuracy difference", quote it off
both runs, not run 1.

Position in the field changed on F and not on H: pydegensac is now **not
separable from poselib-prosac** (-0.0169, CI contains zero) where `master` was
-0.0300 and significant, at 1.27x the leader's cost instead of 3.3x. On H it is
still -0.021 and significant, but at 2.34 ms/pair it is the second-cheapest in
the roster behind cv2-magsac.

### 3b. The original note

The M1 numbers will not transfer -- the two biggest wins were `cov_mat`
(strided passes) and `inlidxs` (a division per correspondence), both of which
are ISA- and compiler-sensitive. `benchmarks/` runs as before; the
equal-correctness arm (`master` + the fixed `lapwrap.c`) is no longer needed
now that the fix is in the branch's history -- compare `08464ca` against
`15445ef` directly.

### 4. Verify output-preservation on Linux -- DONE, gate passes

`scripts/stat_ab.py` over the `pip --target` trees `run_ab.sh` built
(`08464ca` vs `6e72f7f`), 10 golden pairs x 4 metrics x 1000 seeds:

```
min p = 0.03776  (golden_f_reichstag_05461164.../gt_err)  over 40 comparisons
```

The gate is `FAIL_P = 0.01`, so this passes with room. Three of the 40 land
below 0.05 -- `gt_err` and `gt_prec` on one reichstag pair (both 0.03776) and
`inliers` on `golden_h_adam` (0.04279). Two of those three are the same
underlying quantity, since `gt_err` and `gt_prec` are both functions of the
inlier mask, so it is really two independent flags where 40 comparisons at
0.05 predict two.

Worth saying plainly: min p here (0.038) is lower than the 0.108 the docstring
records for `4430bc7`. That is expected rather than worrying -- `4430bc7` was
one reseed removal, this is the entire branch against `master`, including a
different RNG stream and every numerics change in the optimisation pass. The
question these KS tests answer is whether the output *distributions* differ,
not whether the outputs are identical; they are not identical and were never
meant to be (the F inlier checksum moved 11576 -> 11578 across the series).

Note `scipy` is not in `benchmarks/.ab/env` by default -- `pip install scipy`
into whichever interpreter you point at this.

### 4b. The original note

`scripts/stat_ab.py` is now a committed tool rather than an ad hoc script: two
`pip --target` builds in, KS tests out, over inlier count, GT error, GT
precision **and** `model_err` (the returned model scored against ground truth).
Use it the way `f349a6c` was verified. Two things it learned the hard way, both
documented in the file: the first three metrics are functions of the inlier
mask alone and are blind to a model that moves without the mask moving, and the
samples are heavily tied, so values are quantised to nine significant digits
before the KS test -- without that, a 1e-13 shift reads as KS 0.86.

### Traps worth inheriting

- **`PYDEGENSAC_CMAKE_ARGS` persists in the CMake cache.** It silently produced
  a QR-vs-QR comparison here before being caught. `rm -rf build/temp.*` between
  arms -- this applies directly to `benchmarks/toolchain/`, which varies flags
  per arm.
- **The `-O3 -ftree-vectorize -funroll-loops` line in CMakeLists is inert.** It
  sets `CMAKE_CXX_FLAGS` under `CMAKE_COMPILER_IS_GNUCXX`, so it reaches
  `bindings.cpp` on GCC only and has never touched the C core. Measuring those
  flags on the C core showed nothing on M1, so it was left alone -- but do not
  assume the C core has ever been unrolled or vectorised by request.
- **The QR null-space path (`USE_QR`) is repaired but still off.** It was
  unusable: uninitialised `info`, and a `ptrdiff_t[9]` pivot array that LAPACK
  fills with `int32`s, which took the process down with SIGBUS. Fixed in
  `d388f6d` and measured -- 16% slower than the LU path, accuracy tied -- so
  the default is unchanged. Do not spend time re-testing it.

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
