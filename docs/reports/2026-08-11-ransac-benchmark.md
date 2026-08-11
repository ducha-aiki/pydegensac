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

In context: on F, pydegensac's accuracy is **statistically tied with poselib's**
— but it needs 2.9x the time to get there. On H it is measurably *behind* the
leading three while costing 5x what cv2's MAGSAC does. The speed-up narrows
pydegensac's cost gap; it does not put it on the efficient frontier.

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
RANSAC tutorial), RootSIFT-8k mutual-NN pools, ~1500-1900 correspondences per
pair after ratio filtering. mAA over 1-10 deg of `max(R_err, t_err)`. Each
method at its own tuned (px, SNN ratio); iteration budget swept 125 -> 50k.

![F time-mAA](../../benchmarks/results/time_maa_f.png)

| method | best mAA | at budget | mean ms/pair | d mAA vs poselib (paired) |
|---|---|---|---|---|
| poselib | 0.4392 | 50000 | 39.1 | leader |
| poselib-prosac | 0.4363 | 50000 | 106.4 | -0.0027 (-0.0275, +0.0208) |
| pydegensac (base) | 0.4307 | 50000 | 137.1 | -0.0083 (-0.0335, +0.0185) |
| pydegensac (branch) | 0.4297 | 50000 | 111.9 | -0.0097 (-0.0333, +0.0160) |
| cv2-magsac | 0.3762 | 50000 | 28.4 | -0.0633 (-0.0892, -0.0370) `*` |
| cv2-ransac | 0.3403 | 50000 | 282.6 | -0.0987 (-0.1268, -0.0713) `*` |

- **Three-way tie at the top on accuracy**: poselib, poselib-prosac and
  pydegensac are not separable. The two cv2 estimators are, decisively.
- **The separation is in cost.** poselib buys that same accuracy in 39 ms
  against pydegensac's 112 ms. That is the finding that matters, and it is far
  outside timing scatter.
- cv2-ransac is both the least accurate and by far the most expensive here
  (283 ms) — it has no early-termination advantage at these inlier ratios and
  simply runs its budget out.

base vs branch at equal budget:

| budget | mAA base | mAA branch | d mAA (95% CI, paired) | ms base | ms branch | speedup |
|---|---|---|---|---|---|---|
| 125 | 0.3707 | 0.3790 | +0.0081 (-0.0137, +0.0287) | 5.57 | 5.66 | 0.99x |
| 500 | 0.3967 | 0.4055 | +0.0086 (-0.0115, +0.0293) | 7.85 | 7.72 | 1.02x |
| 2500 | 0.4038 | 0.4225 | +0.0183 (-0.0017, +0.0385) | 15.03 | 13.55 | 1.11x |
| 10000 | 0.4168 | 0.4258 | +0.0089 (-0.0100, +0.0290) | 36.54 | 30.66 | 1.19x |
| 25000 | 0.4195 | 0.4287 | +0.0088 (-0.0100, +0.0285) | 75.97 | 62.03 | 1.22x |
| 50000 | 0.4307 | 0.4297 | -0.0014 (-0.0203, +0.0170) | 137.09 | 111.93 | 1.22x |
| **total** | | | | **317.1** | **267.0** | **1.19x** |

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
on macOS. Two things make those an upper bound:

1. **Platform.** The optimisation replaced libc `srandom()`/`random()`. On
   macOS each `srandom()` re-derives the 31-word TYPE_3 state and discards 310
   warm-up draws *behind a lock* (~3.5 us/iteration, profiled at 77% of H
   runtime). glibc runs the same algorithm without the lock, so the cost being
   removed is much smaller here.
2. **Problem size.** The golden pairs are small. These pools carry ~1500-1900
   correspondences, so per-iteration model fitting and error evaluation
   dominate and a fixed per-iteration saving is a smaller share of the whole.

Both show up directly in the data: the speedup is ~1.0x at 125 iterations,
where per-call overheads dominate, and reaches 1.22x only at 50k, where the
per-iteration term is essentially the whole runtime.

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

so ~25% of the gap is the vendored LAPACK, and the rest is still unexplained.
Ruled out by measurement, not assumption: **compiler version** (master built
locally with gcc 10.4 runs at 21.7 ms, indistinguishable from gcc 13.3's
21.5 — so manylinux's gcc 10.2 is not the cause), **assertions** (`NDEBUG` is
set in both, no `__assert_fail`), **hardening flags** (identical `__stack_chk`
usage, no fortify symbols), and **codegen** (both binaries are essentially
all-scalar with near-identical instruction mixes: 59,915 vs 60,779).

Worth pursuing separately from this PR: a newer manylinux image, or linking
OpenBLAS instead of reference LAPACK, looks like free double-digit percent for
every Linux `pip install` user. The residual needs its own investigation.

**Unrelated hazard found on the way**: `pydegensac==0.1.2` silently returns
*every* correspondence as an inlier under numpy 2.x — 300/300 on a synthetic
set with 150 planted outliers, where the same wheel under numpy 1.26 correctly
returns 150. This is the old-pybind11 problem that `0.2` was yanked for, but
**0.1.2 is not yanked**, and it is the version pinned by anything installed
before 2026. It fails silently, which is the worst way to fail.

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
