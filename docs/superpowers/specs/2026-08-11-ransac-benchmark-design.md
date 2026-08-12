# RANSAC speed/accuracy benchmark — design (2026-08-11)

## Purpose

The `speed-up3` branch claims 4.2x (H) / 2.2x (F) on the golden pairs. Those
pairs are five hand-picked images per estimator — enough to gate correctness,
useless as a claim about real workloads. This benchmark answers two questions
on public data with a published protocol:

1. **Does the speed-up survive at scale, and at what accuracy cost?**
   Same-quality-lower-time is the claim; the benchmark has to show quality is
   actually unchanged, not just assert it from the 1000-seed A/B.
2. **Where does pydegensac sit against the current field?** A 2x win is only
   interesting relative to `cv2.USAC_MAGSAC` and poselib, which are what
   people would otherwise use.

Deliverable: a self-contained `benchmarks/` tree anyone can run, plus a report,
a README section and PR text.

## Non-goals

Not a new tuning study, not a claim about pydegensac's *algorithm* being
competitive, and not a general-purpose harness. Thresholds come from the
imc2021-simple study; this repo only reproduces and extends the roster.

## Protocol — F (fundamental)

Port of the `ransac_time_maa` sweep in imc2021-simple
(`docs/reports/2026-07-27-baseline.md`, "RANSAC time/mAA curves").

- **Data**: `RANSAC-Tutorial-Data-ValOnly.tar` → `val/st_peters_square`
  (CVPR-2020 RANSAC tutorial = IMC-2020 PhotoTourism val). 4950 pairs,
  RootSIFT-8k mutual-NN pools, ~2k correspondences/pair; `match_conf` is the
  raw SNN ratio (lower = better), unfiltered (observed range 0.24–1.00).
  Fetched with an HTTP range request (first 600 MB of the 2.2 GB tar covers
  `matches.h5`, `match_conf.h5`, `K1_K2.h5`, `R.h5`, `T.h5`); images are not
  needed — the metric is pose error from GT calibration.
- **Splits**: two disjoint seeded pair subsets — 300 pairs for `tune`,
  600 pairs for evaluation (600 matches the imc21 run).
- **Per method**: filter correspondences at the method's tuned `ratio_th`,
  estimate at its tuned `px_th`, sweep
  `max_iters ∈ {125, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000}`.
- **Metric**: `F` + GT `K1,K2` → `E` → `cv2.recoverPose` against GT relative
  pose; error = max(R_err, t_err) in degrees; mAA = mean over thresholds
  1..10° (imc21 `IMC_DEGREES_THS`, identical to the tutorial's `calc_mAA_FE`).
- **Timing**: `perf_counter` around the estimator call only; report mean and
  median ms/pair. `OMP_NUM_THREADS=1`, `cv2.setNumThreads(1)`.
- **Thresholds**: the imc21 tuned config verbatim —

  | method | ratio_th | px_th |
  |---|---|---|
  | pydegensac | 0.65 | 0.25 |
  | cv2-magsac | 0.65 | 0.25 |
  | poselib, poselib-prosac | 0.70 | 0.50 |

  tuned on **reichstag**, so evaluating on st_peters_square is cross-scene —
  no tuning-on-test. `cv2-ransac` has no imc21 entry and gets one from the
  in-repo `tune` step on the 300-pair tuning subset, over the same grid imc21
  used (`px ∈ {0.25, 0.5, 0.75, 1, 1.5, 2, 4}` × `ratio ∈ {0.6..0.85}`).
  `tune` can be run for every method to reproduce the imc21 values as a check.

## Protocol — H (homography)

Port of the `time_maa` sweep in ds-sac (`bench/time_maa.py`, `bench/report.py`).

- **Data**: `homography.tar.gz` (798 MB) → `EVD` and `HPatchesSeq`, both
  splits. Images are needed here: the metric integrates reprojection error
  over the jointly visible area, which depends on image size.
- **Ground truth for the test split**: the tutorial ships `Hgt.h5` for `val`
  only ("Test(without GT)"). It is recovered from the source datasets:
  - HPatches: `hpatches-sequences-release.zip` (author's HF mirror,
    `vbalnt/hpatches`). Key `{seq}_1_{n}` → `{seq}/H_1_{n}`. The zip central
    directory is read over HTTP range requests, so only the 580 homography
    files are transferred (~163 KB), not the 1.28 GB archive.
  - EVD: `EVD.zip` (28 MB, cmp.felk.cvut.cz). Key `{name}` → `h/{name}.txt`.

  **Verified**: on the `val` split, where the tutorial does ship GT, the
  recovered matrices are bit-identical to `Hgt.h5` for all 145 HPatches pairs
  and all 7 EVD pairs (max |diff| = 0.0). `setup_data.py` re-runs this check
  every time and refuses to write test GT if it fails. HPatches val and test
  sequence sets are disjoint (29 sequences each, no overlap).
- **Splits**: tune on `val` (145 + 7 pairs), report on `test` (145 + 8) —
  the tutorial's own protocol
  (`tune_hyperparameters_and_create_test_H_submission.py`).
- **Per method**: `tune` sweeps `px ∈ {0.25 … 64}` on val and pins one
  threshold per method **per dataset**; test runs that threshold across the
  budget ladder `max_iters ∈ {10, 25, 100, 400, 1600, 6400, 25000}`.
  Correspondences are used as shipped (already SNN ≈ 0.85 pre-filtered), so
  unlike F there is no ratio knob.

  Two deviations from the reference grids, both forced by what tuning found:
  the grid runs to 64 px rather than ds-sac's 4 (four of five methods were
  still improving at 4 px on HPatches, so the reference cap would have
  favoured whichever method peaked inside it), and thresholds are per dataset
  rather than shared (EVD and HPatchesSeq optima differ by 4-16x, and a shared
  threshold ranked on a pooled or averaged score lets EVD's 7 val pairs decide
  what runs on HPatches' 145).
- **Metric**: the tutorial's `get_visible_part_mean_absolute_reprojection_error`
  — reproject image-1 pixel grid under `H_gt` and `H`, average the Euclidean
  residual over the mask of jointly visible pixels; mAA over 10 log-spaced
  thresholds 1–20 px (`calc_mAA`). EVD and HPatches are scored separately and
  never pooled (7-8 vs 145 pairs).

## Methods

Uniform wrapper `fn(pts1, pts2, th, scores) -> (M | None, bool mask)`:

| name | backend |
|---|---|
| `cv2-ransac` | `cv2.findFundamentalMat(FM_RANSAC)` / `cv2.findHomography(RANSAC)` |
| `cv2-magsac` | `cv2.USAC_MAGSAC` |
| `poselib` | `poselib.estimate_fundamental` / `estimate_homography` |
| `poselib-prosac` | same + `progressive_sampling`, input sorted by score |
| `pydegensac` | `pydegensac.findFundamentalMatrix` / `findHomography` |

`conf = 0.9999` (F) / `0.999` (H) for every method, matching imc21 and ds-sac
respectively. PROSAC ordering: poselib does not sort internally, so points are
sorted best-first by `match_conf` and the returned mask is un-permuted — the
wrapper is copied from imc21's `_estimate_fundamental_poselib`. A method that
returns no model, or fewer than the minimum inliers, scores as a failure
(`inf` error), which mAA counts as a miss.

## A/B against main

Only the pydegensac import differs between the two arms, so:

- `run_ab.sh` creates a git worktree at `master` and installs both versions
  with `pip --target` into separate directories, selected per run with
  `PYTHONPATH`. One interpreter and one set of cv2/poselib wheels serve both
  arms, so nothing but pydegensac varies — and neither arm needs rebuilding to
  re-run.
- The branch arm runs **all** methods; the base arm runs `--methods pydegensac`
  only. Halves the runtime, and the other backends are literally the same
  objects in both arms.
- Each result file's meta line carries the backend versions and the pydegensac
  install path (both builds report 0.3.0, so the path is what distinguishes
  them); `report.py` labels the two arms `pydegensac (base)` /
  `pydegensac (branch)`.

Note this machine had neither CMake nor LAPACK, so `pip install .` could not
build the extension at all; `benchmarks/README.md` documents the conda prefix
that supplies both.

## Layout

```
benchmarks/
  README.md          protocol, download sizes, how to reproduce, caveats
  setup_data.py      all four downloads, test-GT recovery + val self-check
  data.py            iter_pairs_f(scene) / iter_pairs_h(dataset, split)
  methods.py         wrappers, method registry, budget factories
  metrics.py         F pose_error + mAA; H visible-part error + mAA
  tuned_config.json  imc21 thresholds + in-repo tuned entries
  run.py             subcommands: tune | f | h  -> results/*.jsonl
  report.py          jsonl -> markdown tables + time-mAA figure
  run_ab.sh          two-venv main-vs-branch driver
  data/              downloads (gitignored)
  .ab/               worktree + the two pydegensac builds (gitignored)
  results/           tuning JSONs and figures (committed); raw per-pair jsonl
                     is gitignored — the F sweep alone is ~27k rows per arm
```

`benchmarks/` is not part of the installed package and is not imported by it;
`setup.py` packaging is untouched.

## Reporting

`docs/reports/2026-08-11-ransac-benchmark.md`: time-mAA curves for both
estimators, main-vs-branch speedup table at matched mAA, and the standing of
pydegensac against the field. The README gets the headline table only.

The imc21 reference numbers (600 st_peters pairs is a different subset from
their 600 reichstag pairs, so absolute mAA will differ) are quoted in the
report for context, not merged into our tables.

## Expected outcome, stated up front

The branch's win is a constant-factor reduction in per-iteration cost. On the
time-mAA plane that shifts pydegensac's curve **left, not up**. Against the
imc21 numbers (pydegensac 0.427 @ 72.8 ms vs poselib-prosac 0.482 @ 17.7 ms at
50k iters) a 2.2x time reduction closes part of the horizontal gap and none of
the vertical one. The report says so plainly; the value of the benchmark does
not depend on pydegensac winning.

## Risks

- **Cross-scene threshold transfer** (F): imc21 itself found tuning transfer
  can fail across datasets. Mitigated by running `tune` on the disjoint
  st_peters tuning subset for at least one method and reporting whether the
  reichstag-tuned values reproduce.
- **Timing noise**: single process, no concurrency, methods interleaved per
  pair so machine drift hits all methods equally; mean *and* median reported.
- **cv2/poselib version drift**: recorded in the jsonl meta so the numbers are
  attributable.
