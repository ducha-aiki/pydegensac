# RANSAC speed/accuracy benchmark (2026-08-11, branch `speed-up3`)

The `speed-up3` session measured 4.2x (H) / 2.2x (F) on five golden pairs per
estimator, on an M-series mac. This benchmark asks what that is worth on public
data at scale, and where pydegensac sits against the estimators people would
otherwise reach for.

Harness: `benchmarks/` (self-contained — `python setup_data.py`, `./run_ab.sh`).
Protocol, data provenance and threshold sourcing: `benchmarks/README.md`.
Machine: WSL2 Linux, glibc, single-threaded, cv2 4.13.0, poselib 2.0.5;
base = `master@08464ca`, branch = `speed-up3@a8caab8`.

## Headline

**The speed-up is real, costs nothing in accuracy, and is much smaller here
than on the golden pairs: 1.19x (F) and 1.16x (H) in aggregate, rising to
1.22x / 1.39x at the largest budgets.** No accuracy difference between the two
builds is significant at any budget on either problem (paired bootstrap over
pairs, 95% CI, all intervals contain zero).

In context, pydegensac is measurably behind the leader on **both** problems:
-0.023 mAA against poselib-prosac on F at 1.5x its cost, and -0.012 against it
on H at 5x what cv2's MAGSAC costs. The speed-up narrows the cost gap; it does
not put pydegensac on the efficient frontier.

**Separately, and much more importantly for macOS users: on macOS the LAPACK
calls were never compiled in at all** — see the section below. Everything in
this report is Linux, where they are, except the M1 section added on
2026-08-12, which re-runs the whole benchmark on Apple silicon against the
fixed build: 1.35x (F) / 1.18x (H) there, and the bug was costing macOS 0.11
mAA on F and 0.21 on H.

## How to read the numbers

Backends are unseeded and the pair sets are finite, so two runs of the
identical configuration differ. A repeat of the whole branch arm gave
|d mAA| up to 0.019 on F and 0.037 on HPatchesSeq, while timings reproduced to
within 1-5%. **Timing differences of a few percent are real; mAA differences
of a couple of points are not.**

Every comparison below is therefore a **paired bootstrap**: all configurations
are scored on the same pairs, so resampling the pair set jointly cancels the
shared pair difficulty and tests the difference directly. `*` marks a 95% CI
that excludes zero. The marginal CIs in the tables are much wider than the
paired ones and should not be used to compare methods.

## F — fundamental matrix

600 pairs of `st_peters_square` (IMC-2020 PhotoTourism val via the CVPR-2020
RANSAC tutorial), RootSIFT-8k mutual-NN pools of ~2000 putative matches, of
which each method's tuned ratio filter keeps a few hundred (median 272 at
ratio 0.85). mAA over 1-10 deg of `max(R_err, t_err)`. Each
method at its own tuned (px, SNN ratio); iteration budget swept 125 -> 50k.

![F time-mAA](../../benchmarks/results/time_maa_f.png)

| method | best mAA | at budget | mean ms/pair | d mAA vs poselib-prosac (paired) |
|---|---|---|---|---|
| poselib-prosac | **0.4570** | 50000 | 40.6 | leader |
| poselib | 0.4392 | 50000 | 39.7 | -0.0178 (-0.0413, +0.0060) |
| pydegensac (branch) | 0.4340 | 25000 | 61.9 | -0.0230 (-0.0452, -0.0003) `*` |
| pydegensac (base) | 0.4307 | 50000 | 137.1 | -0.0261 (-0.0498, -0.0025) `*` |
| cv2-magsac | 0.3762 | 50000 | 28.9 | -0.0811 (-0.1088, -0.0548) `*` |
| cv2-ransac | 0.3403 | 50000 | 288.3 | -0.1165 (-0.1425, -0.0898) `*` |

- **poselib-prosac leads, and pydegensac is measurably behind it** (-0.023,
  paired CI excludes zero) at 1.5x the cost. poselib without PROSAC is not
  separable from the leader.
- The two cv2 estimators are decisively behind, and cv2-ransac is also by far
  the most expensive here (288 ms) — it has no early-termination advantage at
  these inlier ratios and simply runs its budget out.

base vs branch at equal budget:

| budget | mAA base | mAA branch | d mAA (95% CI, paired) | ms base | ms branch | speedup |
|---|---|---|---|---|---|---|
| 125 | 0.3707 | 0.3728 | +0.0020 (-0.0180, +0.0223) | 5.57 | 5.61 | 0.99x |
| 500 | 0.3967 | 0.4042 | +0.0074 (-0.0133, +0.0285) | 7.85 | 7.62 | 1.03x |
| 2500 | 0.4038 | 0.4258 | +0.0220 (+0.0017, +0.0418) `*` | 15.03 | 13.59 | 1.11x |
| 10000 | 0.4168 | 0.4337 | +0.0169 (-0.0027, +0.0365) | 36.54 | 30.61 | 1.19x |
| 25000 | 0.4195 | 0.4340 | +0.0143 (-0.0048, +0.0333) | 75.97 | 61.89 | 1.23x |
| 50000 | 0.4307 | 0.4277 | -0.0032 (-0.0245, +0.0172) | 137.09 | 111.99 | 1.22x |
| **total** | | | | **317.1** | **266.7** | **1.19x** |

One of nine budgets shows a significant accuracy difference (2500, favouring
the branch); with nine comparisons that is what chance produces, and the sign
is not consistent across the ladder.

**The speedup grows monotonically with the iteration budget** (0.99x -> 1.22x).
That is the signature of removing a *per-iteration fixed cost*, which is what
the RNG work did, and it means the saving is largest exactly where pydegensac
is most expensive.

## H — homography

EVD (8 test pairs) and HPatchesSeq (145), scored separately. mAA over 10
log-spaced thresholds 1-20 px of the mean reprojection error over the jointly
visible area. Thresholds tuned per dataset on `val`, reported on `test`.

![H time-mAA](../../benchmarks/results/time_maa_h.png)

HPatchesSeq:

| method | best mAA | at budget | mean ms/pair | d mAA vs poselib-prosac (paired) |
|---|---|---|---|---|
| poselib-prosac | 0.9297 | 1600 | 7.74 | leader |
| cv2-magsac | 0.9269 | 25000 | 0.91 | -0.0027 (-0.0124, +0.0076) |
| poselib | 0.9262 | 25000 | 5.42 | -0.0034 (-0.0214, +0.0097) |
| pydegensac (branch) | 0.9172 | 25000 | 4.53 | -0.0124 (-0.0186, -0.0062) `*` |
| cv2-ransac | 0.9103 | 25000 | 22.65 | -0.0192 (-0.0297, -0.0090) `*` |
| pydegensac (base) | 0.9103 | 25000 | 6.29 | -0.0194 (-0.0290, -0.0110) `*` |

- **The top three are statistically tied; cv2-magsac wins on cost by a mile**
  — 0.91 ms/pair against poselib-prosac's 7.7 and poselib's 5.4, for a
  difference in mAA that the paired test cannot distinguish from zero.
- **pydegensac is measurably behind all three** (-0.012 vs the leader, CI
  excludes zero) at 5x cv2-magsac's cost. The branch does lift it clear of
  master (-0.012 vs -0.019 against the leader).
- base vs branch: 1.05x at small budgets rising to **1.39x at 25000**,
  1.20x aggregate — the same budget-dependent shape as F. No per-budget
  accuracy difference is significant.
- pydegensac is the only method whose HPatches optimum is a tight threshold
  (4 px); every other method wanted 16-64 px, and pydegensac degrades sharply
  when loosened (val mAA 0.9241 at 4 px -> 0.8331 at 16 px). It has less
  headroom from threshold tuning than the rest of the field.

**EVD is not usable for ranking.** With 8 pairs, one pair is 0.125 mAA and the
marginal CIs span 0.20-0.69. Its aggregate speedup reads 1.51x but per-budget
values swing between 0.59x and 1.80x. It is reported because the reference
evaluation includes it, not because it separates anything.

## Why the golden-pair numbers don't transfer

The session report's 4.2x / 2.2x summed per-pair medians over five golden pairs
on macOS. **The platform is the whole story.** The optimisation replaced libc
`srandom()`/`random()`; on macOS each `srandom()` re-derives the 31-word TYPE_3
state and discards 310 warm-up draws *behind a lock* (~3.5 us/iteration,
profiled at 77% of H runtime), while glibc runs the same algorithm without the
lock. The cost being removed here is roughly an order of magnitude smaller.

It is *not* problem size, which was the other candidate: after each method's
tuned ratio filter the estimators see a few hundred correspondences (median 272
at ratio 0.85, 185 at 0.80), not the ~2000 in the raw pools — comparable to the
golden pairs, so that explanation does not survive contact with the data.

The speedup being ~1.0x at 125 iterations and 1.22x at 50k is the other half of
the evidence: a fixed per-iteration saving is invisible when per-call overheads
dominate and worth the most when the iteration term is the whole runtime.

`benchmarks/rng_cost.c` measures the removed cost directly, and confirms the
mechanism quantitatively. On this machine (x86-64, glibc):

| | before (libc reseed + draws) | after (local draws only) | saved |
|---|---|---|---|
| F, 7 draws/iteration | 418.0 ns | 11.8 ns | 406.1 ns (35x) |
| H, 4 draws/iteration | 372.6 ns | 6.8 ns | 365.9 ns (55x) |

406 ns/iteration predicts **20.3 ms** saved over 50k F iterations; the
benchmark measured **25.2 ms** (137.1 -> 111.9). The residual is the rejection
sampling that draws more than 7 values per iteration, plus the LO loops — the
model accounts for ~80% of the observed saving from first principles.

Run the same binary on an M-series mac to settle platform vs. ISA: if the
per-iteration saving there is several times 406 ns, the macOS libc explains the
gap and ARM has nothing to do with it.

This does not diminish the change — 15-20% free on Linux, more on macOS, at an
unchanged output distribution — but the 4.2x / 2.2x figures should be quoted as
macOS golden-pair numbers, not as a general speedup.

## Did anything get lost on Linux along the way? (published releases)

`benchmarks/run_releases.sh` runs pydegensac alone across published PyPI wheels
and local builds of the same code, all on one pinned numpy-1.26 interpreter.
A wheel differs from a local build in *two* ways at once — source and build
environment — so the `v_0.2.2` tag is also built here as the control.

![releases, H](../../benchmarks/results/releases_h.png)

The three cool curves are PyPI wheels, the three warm ones local builds. On
HPatchesSeq they trace the same accuracy at visibly different cost.

Total mean ms/pair summed over the budget ladder:

| arm | F `st_peters_square` | H `HPatchesSeq` |
|---|---|---|
| pypi 0.1.2 (wheel, 2020) | 315.5 | 26.8 |
| pypi 0.2.1 (wheel) | 319.9 | 26.4 |
| pypi 0.2.2 (wheel, latest) | 338.2 | 27.9 |
| local build, tag `v_0.2.2` | 306.5 | 20.9 |
| local build, `master` | 311.2 | 21.5 |
| local build, this branch | **265.4** | **18.1** |

**No regression across releases.** Every published version is within a few
percent of every other on both problems, `master` matches the `v_0.2.2` tag
built the same way, and all six arms are statistically identical in accuracy
(no paired CI excludes zero). Nothing was lost between 0.1.2 and master.

**But the published Linux wheels are slower than the same source built
locally** — 1.10x on F, 1.30-1.34x on H. That is not a regression; it is a
standing property of how the wheels are built. Verified by interleaved
re-measurement (wheel 27.6 / 28.5 / 27.9 ms against local 21.4 / 21.2 /
20.9 ms), so it is not run ordering.

Partially explained. CI builds Linux wheels in manylinux2014 against
`yum lapack-devel`, and auditwheel vendors the result — a **reference
LAPACK/BLAS 3.4.2 from 2012**, against `libgfortran.so.3`. Forcing a local
build to use exactly those bundled libraries via `LD_PRELOAD` costs 1.8 ms of
the 7.0 ms H gap:

| H, 0.2.2 source | ms |
|---|---|
| local build + OpenBLAS | 20.87 |
| local build + Ubuntu reference LAPACK | 21.31 |
| local build + the wheel's bundled LAPACK 3.4.2 | 22.64 |
| the PyPI wheel itself | 27.91 |

so ~25% of the gap is the vendored LAPACK.

### Where the rest goes: a perf profile

Profiled with `perf` on a driver that calls `findHomography` in a loop and does
nothing else (the benchmark's own metric is 85% of process time and buries the
estimator; threads pinned, since unpinned OpenBLAS spends 20-43% in
`__sched_yield` and swamps everything).

| | wheel 0.2.2 | local 0.2.2 |
|---|---|---|
| wall | 0.71 s | 0.54 s |
| `pydegensac.so` | 64.3% -> **0.456 s** | 57.8% -> **0.312 s** |
| `libc` | 15.1% -> 0.107 s | 18.6% -> 0.100 s |
| LAPACK | 3.5% -> 0.025 s | 3.7% -> 0.020 s |

**The extra time is inside pydegensac's own compiled code** — 1.46x on the same
source — while libc and LAPACK are the same in absolute terms. So this is
codegen, not a library.

What that is *not*, each ruled out by direct experiment rather than argument:

- **LAPACK implementation**: 20.87 ms (OpenBLAS) vs 21.31 (modern reference).
- **Optimisation level**: local `-O2` 0.44 s vs `-O3` 0.43 s.
- **Compiler major version**: conda gcc 10.4 21.7 ms vs gcc 13.3 21.5 ms.
- **Symbol visibility**: the wheel exports 425 dynamic symbols against the
  local build's 121, which looked promising — but rebuilding with
  `-fvisibility=hidden` (6 exported) only bought 5%, 0.51 s vs 0.54.
- **Assertions** (`NDEBUG` set in both), **hardening** (identical
  `__stack_chk`, no fortify symbols), **instruction mix** (59,915 vs 60,779,
  both essentially all-scalar).

The one difference not yet isolated: the wheel records `GCC 10.2.1`, which is
RedHat's *devtoolset* build, whereas the gcc-10 control above was conda's
10.4.0 — a different build of GCC with different defaults and a CentOS 7
baseline. Testing that needs the manylinux image itself, i.e. Docker, which is
not available on this machine. (Incidentally the 0.1.2 wheel records two
compilers, GCC 8.5.0 and 12.1.1.)

So the CI change below is worth the ~25% of the gap that is LAPACK; the
remaining ~75% is a real, reproducible property of the manylinux toolchain that
wants a Docker-equipped machine to finish off.

**Unrelated hazard found on the way**: `pydegensac==0.1.2` silently returns
*every* correspondence as an inlier under numpy 2.x — 300/300 on a synthetic
set with 150 planted outliers, where the same wheel under numpy 1.26 correctly
returns 150. This is the old-pybind11 problem that `0.2` was yanked for, but
**0.1.2 is not yanked**, and it is the version pinned by anything installed
before 2026. It fails silently, which is the worst way to fail.

## macOS never called LAPACK at all

`lapwrap.c` declared and called `dgesvd_` and `dsyev_` only under `#ifdef
_WIN32` or `#ifdef __linux__`. macOS defines neither, so on macOS the
preprocessor removed **every LAPACK call in the file**: `lap_SVD` and `lap_eig`
returned without computing anything, leaving `info = 1` and their outputs
untouched.

That is not cosmetic. Both are on the least-squares path:

- `u2f` / `u2fw` (F least squares) call `lap_eig`, then `singulF`, whose first
  act is `if (lap_SVD(...) != 0) { memcpy(F, identity); return; }` — so every
  least-squares F refit on macOS returned the **identity matrix**.
- `u2h` (H least squares, the `len > 4` branch — i.e. every local-optimisation
  refit) calls `lap_eig` and then copies the first 9 values of an untouched
  covariance matrix as the homography.

In other words, local optimisation — the LO in LO-RANSAC, and the whole point
of DEGENSAC's refinement — has been dead on macOS. RANSAC still returns a model
because the garbage refits score badly and get rejected, so the result falls
back to the best minimal-sample hypothesis.

Reproduced by building `899af7f` with `-U__linux__`, which is exactly what
macOS compiles: the resulting `.so` has zero references to `dgesvd_`/`dsyev_`.
Cost, measured on this benchmark:

| | LAPACK linked (Linux) | LAPACK compiled out (macOS) |
|---|---|---|
| F `st_peters_square`, best mAA | 0.4365 | **0.3272** |
| F, total ms over the ladder | 267.6 | **15.7** |
| H `HPatchesSeq`, best mAA | 0.9117 | **0.7021** |
| H, total ms over the ladder | 18.4 | **7.2** |

macOS loses 0.11 mAA on F and 0.21 on H, and is 17x / 2.6x "faster" precisely
because it is skipping the work. Three consequences:

1. **macOS wheels have been shipping a crippled estimator**, for as long as
   these guards have been there.
2. **The 4.2x / 2.2x golden-pair speed-ups in the session report were measured
   on that build**, where the LAPACK path costs nothing — so they describe a
   configuration no Linux or Windows user has ever run.
3. **The golden baselines were captured on macOS** and therefore encode the
   broken behaviour. Fixing this changes macOS outputs and needs a deliberate
   baseline regeneration.

The fix declares both prototypes unconditionally and drops the `#ifdef`s around
the calls. Verified output-preserving *on Linux* — where the calls were already
compiled in — by comparing seeded outputs across 50 real pairs before and
after: byte-identical. On macOS it is by construction a behaviour change.

The same commit removes the per-call workspace query and `malloc`/`free` from
both wrappers (memoised `lwork` per shape, stack buffer). That was expected to
be a speed-up and **is not** — 18.4 ms before and after on H, 267.6 vs 267.8 on
F. The LAPACK path is simply not hot enough on Linux for it to show. It is kept
because it is tidier and removes an allocation from an inner loop, not because
it buys anything measurable.

## The same benchmark on M1, after the fix (2026-08-12)

Everything above is Linux. The whole benchmark was re-run natively on Apple
silicon once the LAPACK fix landed — first because the bug above was found
there independently, and second because the session's headline speed-ups were
macOS numbers and needed restating against a build that does the work.

Machine: M1 MacBook Air, macOS, single-threaded, conda python 3.13,
cv2 4.13.0, poselib 2.0.5 — same data, same thresholds, same protocol.
Three pydegensac arms, because on macOS the branch now differs from `master`
in *two* ways at once:

| arm | build | LAPACK symbols in the `.so` |
|---|---|---|
| `base` | `master@08464ca` | 0 |
| `basefix` | `master` + `f349a6c`'s `lapwrap.c` | 2 |
| `branch` | `f349a6c` | 2 |

**The Linux `-U__linux__` simulation of the bug was accurate.** Measured
natively, on the real macOS build rather than a simulation of it:

| | Linux, simulated (`-U__linux__`) | M1, native (`master`) |
|---|---|---|
| F best mAA | 0.3272 | 0.3133 |
| H `HPatchesSeq` best mAA | 0.7021 | 0.6924 |

Against `branch` on the same M1 run — F 0.4218, H 0.9062 — the bug was costing
macOS users **0.11 mAA on F and 0.21 on H**.

### The speed-up, at equal correctness

`branch` vs `basefix` — both arms call LAPACK, so the only difference is the
RNG/loop work this branch is actually about:

| | aggregate | at the largest budget |
|---|---|---|
| F | **1.35x** | 1.43x @ 50k |
| H `HPatchesSeq` | **1.18x** | 1.57x @ 25k |
| H `EVD` | 2.67x | 3.16x @ 25k |

Accuracy is unaffected: every paired CI contains zero except F at budget 500
(-0.024), which is one marginal result out of 16 and inside the run-to-run
scatter quantified above.

**These are the macOS numbers, and they replace the 4.2x / 2.2x golden-pair
figures**, which were measured on the LAPACK-dead build. Comparing `branch`
against `master` as-shipped on macOS is not a speed comparison at all — the
branch is 11x *slower* there because it is the arm that computes the local
optimisation.

M1 comes out ahead of Linux's 1.19x / 1.16x, in the direction `rng_cost.c`
predicts: it reports **1093 ns saved per F iteration and 1079 ns per H
iteration** on M1 (63x / 110x), against 406 ns on x86-64/glibc. The macOS
`srandom()` lock is worth ~2.7x the glibc saving, which is the part of the
original macOS-vs-Linux gap that was real; the rest of it was the dead LAPACK.

### Where pydegensac sits on M1

| | best mAA | mean ms/pair | leader |
|---|---|---|---|
| F `st_peters_square` | 0.4218 | 96.1 | poselib-prosac 0.4570 @ 50.7 ms |
| H `HPatchesSeq` | 0.9062 | 5.2 | poselib-prosac 0.9297 @ 7.8 ms |
| H `EVD` | 0.3625 | 4.1 | poselib-prosac 0.4500 @ 0.8 ms |

Same qualitative placement as Linux: competitive on F accuracy at roughly 2x
the leader's cost, behind `cv2.USAC_MAGSAC` on homography on both axes. Curves
in `benchmarks/results/time_maa_{f,h}_m1.png`; the equal-correctness pair is in
`time_maa_{f,h}_m1_equalcorrectness.png`, tables alongside as
`report_{f,h}_m1*.md`.

### Golden baselines

Regenerated on M1 with the fixed build (`scripts/make_golden_data.py`); the
gate is green again in exact mode (33 passed). The qualifying pair set moved:
`05466646-05534141` no longer passes capture sanity with local optimisation
restored, and `06373813-06639257` takes its place. Restoring LO changes which
pairs land on the right side of the sanity threshold, so this is expected
rather than a red flag — but it does mean the macOS baselines before and after
this branch describe different code, not just different numbers.

### Reproducing

`setup_data.py` then `run_ab.sh` covers the two-arm case. The third arm was
built by checking `f349a6c`'s `lapwrap.c` into a `master` worktree:

```bash
git -C <master-worktree> checkout f349a6c -- src/pydegensac/degensac/lapwrap.c
pip install --no-deps --target .ab/pkg-basefix <master-worktree>
PYTHONPATH=.ab/pkg-basefix python run.py f --methods pydegensac \
    --label base:m1-master+lapackfix --out f_m1_basefix.jsonl
```

Once this branch is merged the third arm stops being necessary — `master` will
have the fix.

## The tutorial archives disagree about which way `match_conf` points

PROSAC needs correspondences ordered best-first, so the benchmark has to know
whether a dataset's per-match score is "lower is better" (an SNN ratio) or the
opposite. The CVPR-2020 tutorial data does not answer this consistently, and
assuming the SNN convention throughout produced a poselib-prosac curve on
HPatchesSeq that was worse at low budgets than uniform sampling — the tell that
it was sampling worst-first.

Measured against ground truth (`benchmarks/check_scores.py`, AUC of "low score
ranks GT inliers first", 3 px):

| data | raw AUC | pairs agreeing | stored convention |
|---|---|---|---|
| F `st_peters_square` | 0.777 | 200/200 | lower is better |
| H EVD | 0.768 | 8/8 | lower is better |
| H HPatchesSeq | **0.128** | **1%** | **higher is better** |

`data.py` now normalises orientation at load (`H_SCORE_ASCENDING`), so every
consumer can assume lower-is-better, and `check_scores.py` re-derives the table
from ground truth so the flag cannot silently rot. Only the PROSAC backends
read scores, which is what makes this class of bug dangerous: nothing fails,
the estimator just samples badly.

Effect of the fix on HPatchesSeq: poselib-prosac's tuned threshold moved 32 ->
16 px, its val mAA 0.9352 -> 0.9441, and on test it went from fourth place
(0.9166) to the top of the table (0.9297). No other method reads scores, and
the F and EVD orientations were already correct, so nothing else moved.

## PROSAC on F is working, and per-method ratio tuning hid it

The first version of the F table had poselib-prosac at 106 ms against poselib's
39 ms — PROSAC apparently costing 2.7x for nothing. It is not: the ordering is
correct, and the gap was an artefact of tuning each method's SNN ratio
independently.

Two checks. First, the returned inlier mask really does correspond to the
returned model — if the un-permutation after PROSAC's sort were wrong, it would
not: across pairs the mask agrees with "Sampson distance <= threshold" for
98.5-100% of correspondences, the residue being poselib's own post-hoc
refinement. Second, run both at an *identical* configuration (px 0.5,
ratio 0.85), which removes the pool-size difference:

| budget | poselib mAA | ms | poselib-prosac mAA | ms |
|---|---|---|---|---|
| 125 | 0.3342 | 1.57 | **0.4075** | 2.20 |
| 1000 | 0.4350 | 5.10 | **0.4908** | 6.10 |
| 10000 | 0.4992 | 13.92 | **0.5183** | 14.98 |
| 50000 | 0.5042 | 34.72 | **0.5175** | 34.75 |

(120-pair subset, so not comparable to the headline table's absolute values.)
PROSAC is worth +0.073 mAA at 125 iterations for no extra time at large
budgets — exactly what it is supposed to do.

**The tuning rule was the problem.** The grid argmax gave poselib-prosac ratio
0.90 (mAA 0.4433) over 0.85 (0.4413): a 0.002 difference, far inside
run-to-run scatter, but 0.90 keeps a 1.6x bigger pool (median 474 vs 272
correspondences) and cost 2.7x the runtime. Picking a threshold by unqualified
argmax over a noisy grid buys noise with wall-clock. The config now takes the
cheaper option among statistically tied ones, which affected only this entry —
every other method's argmax was either clear or already the cheaper side.

## Thresholds: the imc21 values do not transfer to this scene

The plan was to borrow the tuned thresholds from the imc2021-simple study.
Tuning locally on a 300-pair subset of `st_peters_square`, disjoint from the
600 evaluation pairs, showed that would have been a mistake: every method's
optimum sits at SNN ratio 0.80-0.90, while imc21's reichstag-tuned values are
0.65-0.70. Running at the borrowed values costs:

| method | imc21 (px, ratio) | mAA there | local best | mAA | gap |
|---|---|---|---|---|---|
| pydegensac | 0.25, 0.65 | 0.2603 | 0.5, 0.80 | 0.4077 | **-0.152** |
| poselib-prosac | 0.5, 0.70 | 0.3407 | 0.5, 0.90 | 0.4433 | -0.101 |
| poselib | 0.5, 0.70 | 0.3500 | 0.5, 0.85 | 0.4477 | -0.098 |
| cv2-magsac | 0.25, 0.65 | 0.2737 | 0.25, 0.85 | 0.3637 | -0.090 |

The penalty is several times the spread between methods, and it is uneven — it
would have cost pydegensac 1.7x what it cost cv2-magsac, i.e. it would have
manufactured a result. All five methods are tuned on this scene instead. This
is the "tuning-transfer risk" imc21 flagged for itself, reproduced on a second
pair of scenes.

Both reference grids also had to be widened, because every method's optimum
started at a boundary: the F ratio grid to 1.0 (from 0.85), and the H threshold
grid to 64 px (ds-sac stopped at 4, the tutorial at 2). Full grids in
`benchmarks/results/tuning_{f,h}.json`; the pre-widening F pass is kept in
`tune_f_stage1.log`.

## Caveats

- One machine, one thread, one run per configuration plus one repeat arm for
  scatter. Absolute timings are not portable; the base/branch comparison is,
  since both arms ran the same pairs through the same process layout.
- The two arms share one interpreter and one set of cv2/poselib wheels; only
  the pydegensac import differs.
- Each method is compared at its own best budget, which flatters every method
  equally but is still a selected maximum.
- `poselib-prosac` is slower than `poselib` on F largely because its tuned
  ratio (0.90 vs 0.85) leaves it a bigger correspondence pool. Each method is
  timed at its own optimum — fair, but not a like-for-like pool comparison.
- H thresholds were tuned with the branch build and reused for both arms
  (same algorithm, so the choice is not build-specific).
