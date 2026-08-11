# RANSAC speed/accuracy benchmark

Measures pydegensac against `cv2.RANSAC`, `cv2.USAC_MAGSAC` and poselib
(with and without PROSAC) on public data, as accuracy-vs-compute curves: each
point is one iteration budget, x is the measured time per pair, y is mAA.

Two problems, each following an established evaluation:

| | data | metric | protocol from |
|---|---|---|---|
| **F** | `st_peters_square`, IMC-2020 PhotoTourism val (CVPR-2020 RANSAC tutorial), RootSIFT-8k mutual-NN pools | pose error via `recoverPose` against GT calibration, mAA over 1-10 deg | imc2021-simple `ransac_time_maa` |
| **H** | EVD + HPatchesSeq (same tutorial) | mean reprojection error over the jointly visible area, mAA over 1-20 px | ds-sac `bench/time_maa.py` |

This directory is standalone: it is not part of the installed package and
nothing in `src/` imports it.

## Setup

```bash
python setup_data.py          # ~1.4 GB, idempotent
```

Downloads, in order: 600 MB of the ValOnly tar via a range request (the F
scene's HDF5 files; images are not needed), the 798 MB homography archive
(images *are* needed — the H metric integrates over the visible image area),
and the ground truth for the H **test** split, which the tutorial does not
ship. That last part reads only the ~580 homography files out of the 1.28 GB
`hpatches-sequences-release.zip` over HTTP range requests, plus the 28 MB
`EVD.zip`, and is self-checked: the same lookup is applied to `val`, where the
tutorial *does* ship ground truth, and must reproduce it exactly (it does —
max |diff| = 0 over all 152 val pairs) before any test ground truth is written.

Everything lands in `data/` (gitignored).

## Run

```bash
./run_ab.sh                   # master vs the current branch, both problems
./run_ab.sh v0.2.2            # or any other base ref
```

`run_ab.sh` builds both pydegensac versions into separate directories with
`pip --target` and selects them per run with `PYTHONPATH`, so one interpreter
and one set of cv2/poselib wheels serve both arms — the only thing that
differs is pydegensac. The base arm runs `--methods pydegensac` only; the
branch arm runs the whole roster.

Individually:

```bash
python run.py f --out f.jsonl              # F sweep
python run.py h --out h.jsonl              # H sweep (test split)
python report.py f results/f.jsonl --plot results/time_maa_f.png
```

`$PYTHON` needs numpy, opencv-python, poselib, h5py, matplotlib, and — to
build the extension — a C++ compiler, CMake and LAPACK. If your interpreter
lacks the latter two:

```bash
conda create -y -p .ab/env -c conda-forge python=3.13 numpy h5py matplotlib \
    lapack cmake make
.ab/env/bin/pip install "opencv-python-headless==4.13.*" poselib
PYTHON=$PWD/.ab/env/bin/python ./run_ab.sh
```

## Thresholds

`tuned_config.json` holds one inlier threshold (and, for F, one SNN ratio
threshold) per method **per dataset** — EVD and HPatchesSeq are different
enough that a shared threshold would let EVD's 7 val pairs decide what runs on
HPatches' 145.

- **F** values are taken verbatim from the imc2021-simple tuning study, which
  grid-searched them on 300 **reichstag** pairs of the same feature pools.
  Evaluating on st_peters_square is therefore cross-scene — the thresholds
  were never fitted to the pairs they are scored on. `cv2-ransac` has no entry
  there and is tuned here instead, on a 300-pair tuning subset **disjoint**
  from the 600 evaluation pairs.
- **H** values are all tuned here, on the `val` split, and scored on `test` —
  the tutorial's own protocol. HPatches' val and test sequence sets are
  disjoint (29 sequences each, no overlap).

The H grid runs to 64 px, well past the reference grids (ds-sac stopped at 4,
the tutorial at 2). Most methods are still improving at 4 px on HPatches, so
the reference caps would have handed an advantage to whichever method happened
to peak inside them.

Re-derive either with `python run.py tune f` / `tune h`; the full grid lands in
`results/tuning_{f,h}.json`, whose `best` block is copied into
`tuned_config.json`. Running `tune f` for every method also cross-checks the
borrowed reichstag values against this scene.

## Where did the speed-up go? (`rng_cost.c`)

The branch's win is a per-iteration fixed cost removed from the RANSAC loops.
`rng_cost.c` times exactly that — libc reseed plus minimal-sample draws against
the local generator — so you can tell whether a platform's speed-up is small
because of its libc or because of problem size:

```bash
cc -O3 -I../src/pydegensac/degensac rng_cost.c \
   ../src/pydegensac/degensac/bsd_random.c -o rng_cost && ./rng_cost
```

On x86-64/glibc it reports ~406 ns saved per F iteration, which predicts 20 ms
over a 50k-iteration budget — against 25 ms actually measured, so the model
accounts for most of the observed saving. macOS should report a substantially
larger number, since its `srandom()` takes a lock.

## Protocol notes

- Timing covers the estimator call only, single-threaded (`OMP_NUM_THREADS=1`,
  `cv2.setNumThreads(1)`), with a warm-up call per method before the timed
  loop. Methods are interleaved per pair so machine drift hits them equally.
  Mean **and** median ms/pair are reported: means are skewed by the tail of
  hard pairs that run to the iteration cap.
- A method that returns no model, or too few inliers, scores as a failure and
  mAA counts it as a miss at every threshold. Failure counts are in the tables.
- `conf` is 0.9999 for F and 0.999 for H, matching the two reference studies.
- PROSAC variants get the per-match SNN ratio as quality. poselib does not sort
  internally, so points go in sorted best-first and the mask is un-permuted.
- EVD (7 val / 8 test pairs) and HPatchesSeq (145 each) are scored separately
  and never pooled — EVD is tiny and much harder, so a pooled number would be
  dominated by HPatches while inheriting EVD's noise.
- cv2 is pinned to 4.13 in the reference environment, the version the borrowed
  F thresholds were tuned under. All backend versions are recorded in each
  result file's meta line.
