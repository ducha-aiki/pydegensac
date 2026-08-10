# pydegensac Big Cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add fixed-seed golden regression tests (F on IMC phototourism pairs, H on EVD pairs), remove all unreachable C code and dead duplicates with bit-equivalent results guaranteed by those tests, then fix the verified audit bugs and wire the test suite into CI.

**Architecture:** Two phases. Phase A builds a safety net: a data-generation script extracts SIFT correspondences (reusing imc-2021-simple's extractor), saves inputs + GT matrices (F, K/R/T, H) + current fixed-seed pydegensac outputs as committed `.npz` files, and a pytest module asserts exact equality against them. Phase B deletes dead code in small caller-verified steps, rebuilding and re-running the golden tests after every step.

**Tech Stack:** Python 3.11 (miniforge env `py311`), CMake + LAPACK build via `setup.py`, pytest, OpenCV SIFT, `wxbs-benchmark` (EVD data), imc-2021-simple (SIFT extractor, IMC val data + calibration).

## Global Constraints

- **Bit-equivalence gate:** after every removal step, rebuild and run `python -m pytest tests -q` — golden tests use `np.testing.assert_array_equal` (exact, no tolerance). A step that changes any golden output must be reverted, not tolerated.
- **Rebuild command** (used throughout; run from repo root, in the `py311` env): `pip install . --force-reinstall --no-deps -q && python -m pytest tests -q`
- **PRESERVE VERBATIM:** the `(IDEA: use two separate heuristics, ...)` comment at `src/pydegensac/degensac/exp_ranF.h:10`, and every `//TODO hybrid scoring down to rFtH?` / `//TODO replace with errs[3]...` comment that lives in *surviving* code. Never delete the IDEA comment.
- **Delete only unreachable code.** No merging/rewriting of live functions in this plan (that risks bit-differences); live-code refactoring is deferred to the speed-up phase.
- **Public API frozen:** `findHomography`, `findFundamentalMatrix`, `findHomography_`, `findFundamentalMatrix_`, `convert_cv2_kpts_to_xyA` signatures and behavior unchanged. Downstream contract (imc-2021-simple): `findFundamentalMatrix(Nx2 float64, px_th=, conf=0.9999, max_iters=)` returning `(F, mask)` with `mask.ravel().astype(bool)` working.
- **Do not touch** `lib/pybind11/` or fix audit bugs here (bindings hardening and the resids-NULL speed-up are separate follow-ups; bug fixes could legitimately change failure-path behavior and would muddy the bit-equivalence signal).
- **Baseline provenance:** golden outputs must be produced by a locally-built current-master module. The `py311` env currently imports `/opt/homebrew/.../site-packages/pydegensac` (PyPI wheel) — Task 1 replaces it with a local build before any capture.
- **Deleting a whole dead function also deletes its header declaration**; headers must stay warning-free.
- Data files are committed under `tests/data/` and must stay small (< 2 MB total); tests run offline with no network access.

### Known reachability map (verified 2026-08-10, re-verify with the grep steps before deleting)

Live entry chain from `bindings.cpp`:
- `exp_ransacHcustomLAF` (exp_ranH.c) → `exp_iterHcustom`, `exp_inHranicustom`, `hMCEscustom`, `HcloseToSingular`
- `exp_ransacFcustomLAF` (exp_ranF.c) → `exp_iterF` (**the plain variant — live!**), `exp_inFranicustom`; `HASH_TABLE_F`/`hash.h` are used inside `exp_iterF` → **hash.h stays in exp_ranF.c**
- `DegUtils.c` (F degeneracy path) → `inHrani` (ranH.c) → `iterH` (ranH.c) → **both stay**

Dead (no callers outside themselves / other dead code):
- `ranF.c` — entire file (315 lines: `iterF`, `inFrani`, `ransacF`, `ransacFsimple`, `prosacF`)
- `ranH2el.c` — entire file (549 lines; only external caller of nothing; itself calls `iterH`)
- `ranH.c` — `ransacH`, `ransacHsimple` (keep `inHrani`, `iterH` + any statics they call)
- `exp_ranH.c` — `exp_iterH`, `exp_inHrani`, `hMCEs` (plain copies of the live `*custom` versions)
- `exp_ranF.c` — `exp_ransacF` (line ~256), `exp_ransacFcustom` (line ~836) + helpers exclusively theirs (`exp_inFrani` if only they call it — verify)
- `matutls` — only `mmul, minv, trnm, mattr, rmmult, svduv` are referenced by degensac code; everything not in their transitive closure is dead

---

### Task 1: Local baseline build and provenance check

**Files:**
- No repo changes; environment work only.

**Interfaces:**
- Produces: a `py311` env whose `import pydegensac` resolves to a build of the current master checkout. All later tasks assume this.

- [ ] **Step 1: Build and install the local checkout**

```bash
cd /Users/oldufo/dev/pydegensac
pip install . --force-reinstall --no-deps -q
```

- [ ] **Step 2: Verify provenance and record the baseline commit**

```bash
python -c "import pydegensac, numpy; print(pydegensac.__file__)"
git rev-parse HEAD
```
Expected: the printed path is a fresh site-packages install written by the command above (mtime = now; confirm with `ls -l` on the printed path). Note the commit hash — it is stored inside every golden file in Task 2/3.

- [ ] **Step 3: Run the existing suite to confirm a working baseline**

Run: `python -m pytest tests -q`
Expected: all tests in `test_import.py` and `test_ransac_smoke.py` PASS.

---

### Task 2: Golden data — H on EVD pairs

**Files:**
- Create: `scripts/make_golden_data.py` (H part; Task 3 extends it)
- Create: `tests/data/golden_h_<pair>.npz` (3-5 files, generated)

**Interfaces:**
- Consumes: `extract_sift(image_path: Path, n_features: int, relax_det_th: bool=False) -> (kp (N,2) float32, desc (N,128) float32)` from `/Users/oldufo/dev/imc-2021-simple/examples/sift_eval.py`.
- Produces: `.npz` files with keys `pts1, pts2` ((M,2) float64), `H_gt` (3,3), `H_seed42, mask_seed42, H_seed2024, mask_seed2024`, `px_th, conf, max_iters` (scalars), `baseline_commit` (bytes). Task 4's test reads exactly these keys.

- [ ] **Step 1: Install the data dependency**

```bash
pip install wxbs-benchmark
python -c "from wxbs_benchmark.dataset import EVDDataset; print('ok')"
```
Expected: `ok`.

- [ ] **Step 2: Write `scripts/make_golden_data.py` with the shared helpers and the H part**

```python
#!/usr/bin/env python3
"""Generate golden regression data for pydegensac.

Inputs are SIFT correspondences (extracted with imc-2021-simple's
extractor), ground truth comes with the datasets (EVD homographies,
IMC calibration), and the "golden" outputs are the current pydegensac
results at fixed seeds.  Regenerate ONLY deliberately: committing new
goldens redefines the bit-equivalence baseline.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import cv2
import numpy as np

REPO = Path(__file__).resolve().parent.parent
DATA_DIR = REPO / "tests" / "data"
IMC_SIMPLE = Path("/Users/oldufo/dev/imc-2021-simple")
sys.path.insert(0, str(IMC_SIMPLE / "examples"))
sys.path.insert(0, str(IMC_SIMPLE))
from sift_eval import extract_sift  # noqa: E402

import pydegensac  # noqa: E402

N_FEATURES = 8000
RATIO_TH = 0.8
SEEDS = (42, 2024)
H_PARAMS = dict(px_th=3.0, conf=0.999, max_iters=20000)
F_PARAMS = dict(px_th=0.75, conf=0.9999, max_iters=10000)

BASELINE_COMMIT = subprocess.check_output(
    ["git", "rev-parse", "HEAD"], cwd=REPO).decode().strip()


def ratio_match(kp_a, desc_a, kp_b, desc_b, ratio_th=RATIO_TH):
    """Lowe-ratio BF matching, mirroring imc-2021-simple match_descriptors."""
    if len(kp_a) < 8 or len(kp_b) < 8:
        return None, None
    bf = cv2.BFMatcher(cv2.NORM_L2)
    raw = bf.knnMatch(desc_a, desc_b, k=2)
    good = [(m, n) for m, n in raw if m.distance < ratio_th * n.distance]
    if len(good) < 8:
        return None, None
    mkp_a = kp_a[[m.queryIdx for m, _ in good]].astype(np.float64)
    mkp_b = kp_b[[m.trainIdx for m, _ in good]].astype(np.float64)
    return mkp_a, mkp_b


def sift_on_array(img, n_features=N_FEATURES):
    """Run extract_sift on an in-memory image via a temp file (reuses the
    exact imc-2021-simple code path rather than reimplementing it)."""
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
        cv2.imwrite(f.name, img)
        return extract_sift(Path(f.name), n_features=n_features)


def run_h(pts1, pts2, seed):
    H, mask = pydegensac.findHomography(pts1, pts2, seed=seed, **H_PARAMS)
    mask = np.asarray(mask, dtype=bool)
    # determinism check: capture must be reproducible before it is golden
    H2, mask2 = pydegensac.findHomography(pts1, pts2, seed=seed, **H_PARAMS)
    assert np.array_equal(H, H2) and np.array_equal(mask, np.asarray(mask2, dtype=bool))
    return H, mask


def h_sanity(H, mask, pts1, pts2, H_gt, th=10.0):
    """At least half of the found inliers must agree with GT H within th px."""
    if mask.sum() < 10:
        return False
    p1 = np.concatenate([pts1[mask], np.ones((mask.sum(), 1))], axis=1)
    proj = (H_gt @ p1.T).T
    proj = proj[:, :2] / proj[:, 2:3]
    err = np.linalg.norm(proj - pts2[mask], axis=1)
    return np.median(err) < th


def make_evd_goldens(max_pairs=5):
    from wxbs_benchmark.dataset import EVDDataset
    dset = EVDDataset(str(REPO / ".cache" / "EVD"), download=True)
    saved = 0
    for pair in dset:
        keys = set(pair.keys())
        h_key = "H" if "H" in keys else ("H_gt" if "H_gt" in keys else None)
        assert h_key is not None, f"no homography key in EVD sample: {keys}"
        name = pair.get("name", f"pair{saved}")
        kp1, d1 = sift_on_array(pair["img1"])
        kp2, d2 = sift_on_array(pair["img2"])
        pts1, pts2 = ratio_match(kp1, d1, kp2, d2)
        if pts1 is None or len(pts1) < 20:
            print(f"skip {name}: too few matches")
            continue
        H_gt = np.asarray(pair[h_key], dtype=np.float64)
        results = {}
        ok = True
        for seed in SEEDS:
            H, mask = run_h(pts1, pts2, seed)
            if not h_sanity(H, mask, pts1, pts2, H_gt):
                ok = False
                break
            results[f"H_seed{seed}"] = H
            results[f"mask_seed{seed}"] = mask
        if not ok:
            print(f"skip {name}: sanity check failed")
            continue
        out = DATA_DIR / f"golden_h_{name}.npz"
        np.savez_compressed(
            out, pts1=pts1, pts2=pts2, H_gt=H_gt,
            baseline_commit=np.bytes_(BASELINE_COMMIT.encode()),
            **{k: np.float64(v) if not isinstance(v, int) else np.int64(v)
               for k, v in H_PARAMS.items()},
            **results)
        print(f"saved {out.name}: {len(pts1)} matches, "
              f"{results['mask_seed42'].sum()} inliers")
        saved += 1
        if saved >= max_pairs:
            break
    assert saved >= 3, f"only {saved} EVD pairs qualified; relax RATIO_TH to 0.9 and retry"


if __name__ == "__main__":
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    make_evd_goldens()
```

Note: if `EVDDataset` samples use different keys than `img1`/`img2` (inspect the first `pair` with `print(pair.keys())` on first failure), adapt the two access lines — the dataset API is the only part not verifiable offline. If SIFT finds too few matches on the extreme EVD pairs even at ratio 0.9, fall back to the pre-extracted EVD from https://cmp.felk.cvut.cz/wbs/ (download `EVD.zip`; `h/` holds GT homographies as text files, `1/` and `2/` the image pairs) and keep the same npz schema.

- [ ] **Step 3: Run it**

Run: `python scripts/make_golden_data.py`
Expected: `saved golden_h_*.npz` lines for 3-5 pairs, each with >= 10 inliers; script exits 0.

- [ ] **Step 4: Inspect sizes and commit**

```bash
du -ch tests/data/golden_h_*.npz
git add scripts/make_golden_data.py tests/data/golden_h_*.npz
git commit -m "add golden H regression data (EVD, fixed seeds)"
```
Expected: total well under 1 MB.

---

### Task 3: Golden data — F on IMC phototourism pairs

**Files:**
- Modify: `scripts/make_golden_data.py` (add F part)
- Create: `tests/data/golden_f_<scene>_<a>_<b>.npz` (3-5 files, generated)

**Interfaces:**
- Consumes: `download_dataset(dataset, split, dest) -> str` and `list_split_scenes` from `imc2021.download`; `load_calibration(scene_dir) -> dict[name -> {'K','R','T'}]` from `imc2021.io`; `pose_error` from `imc2021.metrics`; `extract_sift`, `ratio_match` from Task 2.
- Produces: `.npz` files with keys `pts1, pts2, F_gt, K1, R1, T1, K2, R2, T2`, `F_seed42, mask_seed42, F_seed2024, mask_seed2024`, `px_th, conf, max_iters`, `baseline_commit`. Task 4's test reads exactly these keys.

- [ ] **Step 1: Add GT-F computation and the IMC part to the script**

Append to `scripts/make_golden_data.py`:

```python
def f_from_krt(c1, c2):
    """GT fundamental matrix from two calibrations (world->cam: x = R X + T)."""
    R = c2["R"] @ c1["R"].T
    t = (c2["T"].flatten() - R @ c1["T"].flatten())
    tx = np.array([[0, -t[2], t[1]],
                   [t[2], 0, -t[0]],
                   [-t[1], t[0], 0]])
    E = tx @ R
    F = np.linalg.inv(c2["K"]).T @ E @ np.linalg.inv(c1["K"])
    return F / F[2, 2] if abs(F[2, 2]) > 1e-12 else F


def run_f(pts1, pts2, seed):
    F, mask = pydegensac.findFundamentalMatrix(pts1, pts2, seed=seed, **F_PARAMS)
    mask = np.asarray(mask, dtype=bool)
    F2, mask2 = pydegensac.findFundamentalMatrix(pts1, pts2, seed=seed, **F_PARAMS)
    assert np.array_equal(F, F2) and np.array_equal(mask, np.asarray(mask2, dtype=bool))
    return F, mask


def f_sanity(F, mask, pts1, pts2, c1, c2, max_deg=15.0):
    from imc2021.metrics import pose_error
    if mask.sum() < 30:
        return False
    err_R, err_t = pose_error(F, pts1[mask], pts2[mask], c1, c2)
    return max(err_R, err_t) < max_deg


def make_imc_goldens(max_pairs=5, scene="reichstag"):
    from imc2021.download import download_dataset
    from imc2021.io import load_calibration, find_image
    root = Path(download_dataset("phototourism", "val", dest=str(REPO / ".cache" / "imc")))
    scene_dir = root / scene
    calib = load_calibration(scene_dir)
    names = sorted(calib.keys())
    images_dir = scene_dir / "set_100" / "images"
    feats = {}
    saved = 0
    # walk consecutive-ish pairs until enough qualify
    for i in range(0, len(names) - 1):
        a, b = names[i], names[i + 1]
        for n in (a, b):
            if n not in feats:
                img_path = find_image(n, images_dir)
                feats[n] = extract_sift(img_path, n_features=N_FEATURES,
                                        relax_det_th=True)
        pts1, pts2 = ratio_match(*feats[a], *feats[b])
        if pts1 is None or len(pts1) < 100:
            continue
        F_gt = f_from_krt(calib[a], calib[b])
        results = {}
        ok = True
        for seed in SEEDS:
            F, mask = run_f(pts1, pts2, seed)
            if not f_sanity(F, mask, pts1, pts2, calib[a], calib[b]):
                ok = False
                break
            results[f"F_seed{seed}"] = F
            results[f"mask_seed{seed}"] = mask
        if not ok:
            continue
        out = DATA_DIR / f"golden_f_{scene}_{a}_{b}.npz"
        np.savez_compressed(
            out, pts1=pts1, pts2=pts2, F_gt=F_gt,
            K1=calib[a]["K"], R1=calib[a]["R"], T1=calib[a]["T"],
            K2=calib[b]["K"], R2=calib[b]["R"], T2=calib[b]["T"],
            baseline_commit=np.bytes_(BASELINE_COMMIT.encode()),
            **{k: np.float64(v) if not isinstance(v, int) else np.int64(v)
               for k, v in F_PARAMS.items()},
            **results)
        print(f"saved {out.name}: {len(pts1)} matches, "
              f"{results['mask_seed42'].sum()} inliers")
        saved += 1
        if saved >= max_pairs:
            break
    assert saved >= 3, f"only {saved} IMC pairs qualified; widen the pair walk or try another scene"
```

and change the `__main__` block to:

```python
if __name__ == "__main__":
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    make_evd_goldens()
    make_imc_goldens()
```

Adapt the exact `find_image`/`pose_error` signatures to what `imc2021` actually exports (check `imc2021/io.py` and `imc2021/metrics.py:36` — `pose_error(F, kp_a, kp_b, calib_a, calib_b)` returns angular errors; if it returns a single dict/tuple shape differing from `(err_R, err_t)`, unpack accordingly). `download_dataset("phototourism", "val")` downloads ~a few GB once into `.cache/imc` — reuse the existing `/Users/oldufo/dev/imc-2021-simple/data` directory instead if the user already has it (pass `dest=` accordingly).

- [ ] **Step 2: Run it**

Run: `python scripts/make_golden_data.py`
Expected: 3-5 `golden_f_reichstag_*.npz` saved (plus the H files regenerated identically — verify `git diff --stat` shows no change to the H files; if it does, the capture is nondeterministic and MUST be debugged before proceeding).

- [ ] **Step 3: Add `.cache` to .gitignore and commit**

```bash
echo ".cache/" >> .gitignore
git add .gitignore scripts/make_golden_data.py tests/data/golden_f_*.npz
git commit -m "add golden F regression data (IMC reichstag, fixed seeds)"
```

---

### Task 4: Fixed-seed regression tests

**Files:**
- Create: `tests/test_golden_regression.py`
- Test: `tests/test_golden_regression.py` (self)

**Interfaces:**
- Consumes: the `.npz` schemas from Tasks 2-3.
- Produces: the bit-equivalence gate every Phase-B step runs.

- [ ] **Step 1: Write the test file**

```python
import os
from pathlib import Path

import numpy as np
import pytest

import pydegensac

DATA_DIR = Path(__file__).parent / "data"
H_FILES = sorted(DATA_DIR.glob("golden_h_*.npz"))
F_FILES = sorted(DATA_DIR.glob("golden_f_*.npz"))
SEEDS = (42, 2024)

# Bit-exact comparison only holds on the platform/toolchain that captured the
# baselines (rand()/random() and FP codegen differ across libc's).  Locally
# this is the refactoring gate (exact); CI sets PYDEGENSAC_GOLDEN_EXACT=0 and
# gets sanity + within-run determinism instead.
EXACT = os.environ.get("PYDEGENSAC_GOLDEN_EXACT", "1") == "1"


def _assert_matches_golden(model, mask, d, model_key, mask_key, rerun):
    mask = np.asarray(mask, dtype=bool)
    if EXACT:
        np.testing.assert_array_equal(mask, d[mask_key])
        np.testing.assert_array_equal(model, d[model_key])
    else:
        # cross-platform mode: determinism within this binary + rough agreement
        model2, mask2 = rerun()
        np.testing.assert_array_equal(model, model2)
        np.testing.assert_array_equal(mask, np.asarray(mask2, dtype=bool))
        assert mask.sum() >= max(10, 0.5 * d[mask_key].sum())


def test_golden_data_exists():
    assert len(H_FILES) >= 3, "golden H data missing - run scripts/make_golden_data.py"
    assert len(F_FILES) >= 3, "golden F data missing - run scripts/make_golden_data.py"


@pytest.mark.parametrize("path", H_FILES, ids=lambda p: p.stem)
@pytest.mark.parametrize("seed", SEEDS)
def test_homography_bit_equivalent(path, seed):
    d = np.load(path)
    call = lambda: pydegensac.findHomography(
        d["pts1"], d["pts2"],
        px_th=float(d["px_th"]), conf=float(d["conf"]),
        max_iters=int(d["max_iters"]), seed=seed)
    H, mask = call()
    _assert_matches_golden(H, mask, d, f"H_seed{seed}", f"mask_seed{seed}", call)


@pytest.mark.parametrize("path", F_FILES, ids=lambda p: p.stem)
@pytest.mark.parametrize("seed", SEEDS)
def test_fundamental_bit_equivalent(path, seed):
    d = np.load(path)
    call = lambda: pydegensac.findFundamentalMatrix(
        d["pts1"], d["pts2"],
        px_th=float(d["px_th"]), conf=float(d["conf"]),
        max_iters=int(d["max_iters"]), seed=seed)
    F, mask = call()
    _assert_matches_golden(F, mask, d, f"F_seed{seed}", f"mask_seed{seed}", call)


@pytest.mark.parametrize("path", H_FILES[:1] + F_FILES[:1], ids=lambda p: p.stem)
def test_gt_sanity(path):
    """Loose GT check so regenerated baselines can't silently go bad."""
    d = np.load(path)
    mask = d["mask_seed42"]
    assert mask.sum() >= 10
    if "H_gt" in d:
        p1 = np.concatenate([d["pts1"][mask], np.ones((mask.sum(), 1))], axis=1)
        proj = (d["H_gt"] @ p1.T).T
        err = np.linalg.norm(proj[:, :2] / proj[:, 2:3] - d["pts2"][mask], axis=1)
        assert np.median(err) < 10.0
    else:
        p1 = np.concatenate([d["pts1"][mask], np.ones((mask.sum(), 1))], axis=1)
        p2 = np.concatenate([d["pts2"][mask], np.ones((mask.sum(), 1))], axis=1)
        # symmetric epipolar sanity against GT F
        Fx1 = (d["F_gt"] @ p1.T).T
        num = np.abs(np.sum(p2 * Fx1, axis=1))
        den = np.sqrt(Fx1[:, 0] ** 2 + Fx1[:, 1] ** 2)
        assert np.median(num / den) < 5.0
```

- [ ] **Step 2: Run the new tests**

Run: `python -m pytest tests/test_golden_regression.py -q`
Expected: all PASS (they test the same build that produced the goldens).

- [ ] **Step 3: Full-suite run and commit**

```bash
python -m pytest tests -q
git add tests/test_golden_regression.py
git commit -m "add fixed-seed golden regression tests"
```

---

### Task 5: Delete whole dead files (ranF.c, ranH2el.c, stale build detritus)

**Files:**
- Delete: `src/pydegensac/degensac/ranF.c`, `ranF.h`, `ranH2el.c`, `ranH2el.h`, `src/pydegensac/degensac/Makefile`, `src/pydegensac/matutls/Makefile`, `src/pydegensac/matutls/scons.what`, `src/pydegensac/matutls/matutls` (stray file), `src/pydegensac/matutls/solv.s`, `src/pydegensac/matutls/svd2.c` (not in any build)
- Modify: `CMakeLists.txt:31-44` (remove `ranF.c`, `ranH2el.c` from `degensac_srcs`)

**Interfaces:**
- Consumes: golden gate from Task 4.
- Produces: a build without the two dead translation units.

- [ ] **Step 1: Re-verify the files are unreferenced**

```bash
cd src/pydegensac/degensac
grep -rn "ranF\.h\|ranH2el\.h" *.c *.h ../bindings.cpp | grep -v "^ranF\.\|^ranH2el\.\|exp_ranF"
```
Expected: no output (nothing live includes them). If `exp_ranF.c` includes `ranF.h`, keep `ranF.h` but still delete `ranF.c` — declarations cost nothing, undefined-but-unused symbols link fine in a static lib only if never referenced.

- [ ] **Step 2: Delete and update CMake**

```bash
git rm src/pydegensac/degensac/ranF.c src/pydegensac/degensac/ranF.h \
       src/pydegensac/degensac/ranH2el.c src/pydegensac/degensac/ranH2el.h \
       src/pydegensac/degensac/Makefile src/pydegensac/matutls/Makefile \
       src/pydegensac/matutls/scons.what src/pydegensac/matutls/solv.s \
       src/pydegensac/matutls/svd2.c
git rm src/pydegensac/matutls/matutls 2>/dev/null || true
```
Then edit `CMakeLists.txt` `degensac_srcs`: remove the `ranF.c` and `ranH2el.c` lines.

- [ ] **Step 3: Rebuild + bit-equivalence gate**

Run: `pip install . --force-reinstall --no-deps -q && python -m pytest tests -q`
Expected: builds clean, all tests PASS.

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "remove unreachable ranF.c, ranH2el.c and stale build files"
```

---

### Task 6: Trim ranH.c to the degeneracy-check subset

**Files:**
- Modify: `src/pydegensac/degensac/ranH.c` (keep `inHrani` + `iterH` + statics they call; delete `ransacH`, `ransacHsimple`)
- Modify: `src/pydegensac/degensac/ranH.h` (drop deleted declarations)

**Interfaces:**
- Consumes: `DegUtils.c` calls `inHrani`; `inHrani` calls `iterH` (verified).
- Produces: `ranH.c` exporting only `iterH`, `inHrani`.

- [ ] **Step 1: Verify the delete list has no live callers**

```bash
cd src/pydegensac/degensac
grep -wn "ransacH\|ransacHsimple" *.c *.h ../bindings.cpp
```
Expected: hits only inside `ranH.c`/`ranH.h`. Also list every static/helper `ransacH` uses that `iterH`/`inHrani` do not — delete those too, keep shared ones.

- [ ] **Step 2: Delete `ransacH` and `ransacHsimple` bodies from ranH.c and their declarations from ranH.h**

Keep the file header comment block intact. Do not touch `iterH`/`inHrani` bodies in any way (not even whitespace).

- [ ] **Step 3: Rebuild + gate**

Run: `pip install . --force-reinstall --no-deps -q && python -m pytest tests -q`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add -A && git commit -m "trim ranH.c to the inHrani/iterH subset used by the degeneracy check"
```

---

### Task 7: Delete dead duplicate estimators in exp_ranH.c and exp_ranF.c

**Files:**
- Modify: `src/pydegensac/degensac/exp_ranH.c` (delete `exp_iterH`, `exp_inHrani`, `hMCEs`)
- Modify: `src/pydegensac/degensac/exp_ranH.h` (drop their declarations)
- Modify: `src/pydegensac/degensac/exp_ranF.c` (delete `exp_ransacF`, `exp_ransacFcustom` + exclusively-theirs helpers)
- Modify: `src/pydegensac/degensac/exp_ranF.h` (drop their declarations — **the IDEA comment at line 10 stays**)

**Interfaces:**
- Consumes: live chains from the reachability map (exp_ransacHcustomLAF → `*custom` helpers; exp_ransacFcustomLAF → `exp_iterF` + `exp_inFranicustom`).
- Produces: one estimator implementation per model type.

- [ ] **Step 1: Caller-verify each candidate in exp_ranH.c**

```bash
cd src/pydegensac/degensac
for fn in exp_iterH exp_inHrani hMCEs; do
  echo "== $fn"; grep -wn "$fn" *.c *.h ../bindings.cpp
done
```
Expected: `exp_iterH` referenced only in `exp_ranH.c:54` (def), `exp_ranH.c:228` (inside dead `exp_inHrani`), and `exp_ranH.h`; `exp_inHrani` only def + header; `hMCEs` (exact word, not `hMCEscustom`) only def + header + possibly calls from the two dead functions. Any hit inside `exp_ransacHcustomLAF`, `DegUtils.c`, or `bindings.cpp` → that function is LIVE, remove it from the delete list.

- [ ] **Step 2: Delete the exp_ranH.c dead trio + their header declarations; rebuild + gate**

Run: `pip install . --force-reinstall --no-deps -q && python -m pytest tests -q`
Expected: PASS. Commit: `git add -A && git commit -m "remove dead plain-variant estimators from exp_ranH.c"`

- [ ] **Step 3: Caller-verify the exp_ranF.c candidates**

```bash
for fn in exp_ransacF exp_ransacFcustom exp_inFrani; do
  echo "== $fn"; grep -wn "$fn" *.c *.h ../bindings.cpp
done
```
Expected: no references from `exp_ransacFcustomLAF`'s body (lines ~1281+), `bindings.cpp`, or `DegUtils.c`. **Known trap:** `exp_iterF` IS live (called by `exp_ransacFcustomLAF`) — do not delete it; `HASH_TABLE_F`, `htInit`/`htClear`/`htContains` uses inside surviving functions keep `hash.h` and `hash.c`. If deleting `exp_ransacF` removes the only `htInit(&HASH_TABLE_F)` call: static storage is zero-initialized, so behavior of the live path is unchanged — the golden gate confirms.

- [ ] **Step 4: Delete the dead F estimators + exclusively-theirs helpers + header declarations, preserving exp_ranF.h:10 IDEA comment; rebuild + gate**

Run: `pip install . --force-reinstall --no-deps -q && python -m pytest tests -q`
Expected: PASS, and `grep -c "IDEA" src/pydegensac/degensac/exp_ranF.h` prints `1`.

- [ ] **Step 5: Commit**

```bash
git add -A && git commit -m "remove dead exp_ransacF/exp_ransacFcustom duplicates"
```

---

### Task 8: Trim matutls to the used subset

**Files:**
- Modify: `src/pydegensac/matutls/CMakeLists.txt` (shrink `matutls_srcs`)
- Delete: every `src/pydegensac/matutls/*.c` not in the computed keep-set

**Interfaces:**
- Consumes: entry points `mmul, minv, trnm, mattr, rmmult, svduv` (referenced from degensac C) + whatever LAPACK wrappers `lapwrap.c` needs.
- Produces: minimal `matutls` static lib.

- [ ] **Step 1: Compute the transitive keep-set from the built objects**

Build once (`pip install . -q --force-reinstall --no-deps`), find the CMake build dir under `build/temp.*`, then:

```python
# scripts/_matutls_closure.py  (throwaway; do not commit)
import subprocess, sys, glob, os
from collections import defaultdict

objdir = sys.argv[1]  # .../build/temp.*/src/pydegensac/matutls/CMakeFiles/matutls.dir
objs = glob.glob(os.path.join(objdir, "*.o"))
defined, undefined = {}, defaultdict(set)
for o in objs:
    out = subprocess.check_output(["nm", "-g", o]).decode()
    for line in out.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[-2] in ("T", "D", "S"):
            defined[parts[-1]] = o
        elif parts and parts[0] == "U" or (len(parts) >= 2 and parts[-2] == "U"):
            undefined[o].add(parts[-1])

roots = {"_" + n for n in ["mmul", "minv", "trnm", "mattr", "rmmult", "svduv"]}
keep, work = set(), [defined[s] for s in roots if s in defined]
while work:
    o = work.pop()
    if o in keep:
        continue
    keep.add(o)
    for sym in undefined[o]:
        if sym in defined and defined[sym] not in keep:
            work.append(defined[sym])
print("KEEP:", sorted(os.path.basename(o).replace(".c.o", ".c") for o in keep))
```

Run: `python scripts/_matutls_closure.py <objdir>`
Expected: a keep list of roughly 10-20 files (svduv pulls in the qr/householder chain: expect `qrbdv.c`, `ldumat.c`, `ldvmat.c`, `housev.c` or similar). Always keep `matconsts.c` and `ccmath_globals.c` (global definitions).

- [ ] **Step 2: Rewrite `matutls_srcs` in `src/pydegensac/matutls/CMakeLists.txt` to exactly the keep list + `matconsts.c` + `ccmath_globals.c`; `git rm` all other `.c` files**

Keep all `.h` files (`matutl.h`, `ccmath.h`, `complex.h` — headers may cross-reference).

- [ ] **Step 3: Rebuild + gate**

Run: `pip install . --force-reinstall --no-deps -q && python -m pytest tests -q`
Expected: link succeeds (a missed transitive dep fails loudly at link time — add it back and re-run), all tests PASS.

- [ ] **Step 4: Run the example smoke test end to end**

```bash
cd examples && python simple-example.py && cd ..
```
Expected: homography + fundamental output like the README (requires cv2 — already in env).

- [ ] **Step 5: Commit**

```bash
git add -A && git commit -m "trim matutls to the transitively used subset"
```

---

### Task 9: Final sweep, downstream smoke test, docs

**Files:**
- Modify: `src/pydegensac/degensac/*.c/h` (remove now-dangling includes/decls only — verified by compiler warnings)
- Modify: `CLAUDE.md` (architecture section: dead files are gone; update the file examples)

**Interfaces:**
- Consumes: everything above.
- Produces: clean build, verified downstream compatibility.

- [ ] **Step 1: Warning sweep**

```bash
pip install . --force-reinstall --no-deps -q 2>&1 | grep -i "warning.*unused" | sort | uniq -c | head -30
```
Expected: no unused-function/unused-variable warnings pointing at leftovers of the deletion (pre-existing legacy warnings in live code stay — do not chase them).

- [ ] **Step 2: Full gate one last time, twice (determinism)**

Run: `python -m pytest tests -q && python -m pytest tests/test_golden_regression.py -q`
Expected: PASS both runs.

- [ ] **Step 3: Downstream smoke test against the new build**

```bash
cd /Users/oldufo/dev/imc-2021-simple && python -m pytest tests/test_ransac.py -q -k pydegensac; cd /Users/oldufo/dev/pydegensac
```
Expected: PASS (uses the same `py311` env, which now has the trimmed local build).

- [ ] **Step 4: Update CLAUDE.md architecture bullets**

In the "Architecture" section, layer 1: replace the file-example list (`ranH.c`, `ranF.c`, ...) with the surviving set (`exp_ranH.c`, `exp_ranF.c`, `ranH.c` (degeneracy-support subset), `DegUtils.c`, `Ftools.c`, `Htools.c`, `rtools.c`, `utools.c`, `hash.c`, `lapwrap.c`) and add one line: "Golden fixed-seed regression tests in `tests/test_golden_regression.py` guarantee bit-equivalence; regenerate baselines only deliberately via `scripts/make_golden_data.py`."

- [ ] **Step 5: Commit**

```bash
git add -A && git commit -m "final cleanup sweep, docs, downstream smoke-tested"
```

---

### Task 10: Fix Python-layer bugs (utils.py)

**Files:**
- Modify: `src/pydegensac/utils.py`
- Test: `tests/test_input_validation.py` (created here, extended in Task 11)

**Interfaces:**
- Consumes: public API from `utils.py`.
- Produces: failure path returns `np.ndarray` bool mask; typed exceptions instead of bare `except:`.

- [ ] **Step 1: Write failing tests**

```python
# tests/test_input_validation.py
import numpy as np
import pytest

import pydegensac
from pydegensac import utils


def test_failure_path_mask_is_bool_ndarray(monkeypatch):
    n = 8
    pts = np.random.RandomState(0).rand(n, 2)
    # force the no-model branch by stubbing the raw binding
    monkeypatch.setattr(utils.pydegensac, "findHomography_",
                        lambda *a, **k: (np.zeros((3, 3)), np.zeros(n)))
    H, mask = pydegensac.findHomography(pts, pts + 1)
    assert isinstance(mask, np.ndarray) and mask.dtype == np.bool_
    assert mask.shape == (n,) and not mask.any()
    monkeypatch.setattr(utils.pydegensac, "findFundamentalMatrix_",
                        lambda *a, **k: (np.zeros((3, 3)), np.zeros(n)))
    F, maskf = pydegensac.findFundamentalMatrix(pts, pts + 1)
    assert isinstance(maskf, np.ndarray) and maskf.dtype == np.bool_


def test_bad_error_type_raises_valueerror():
    pts = np.random.RandomState(0).rand(8, 2)
    with pytest.raises(ValueError, match="Error type"):
        pydegensac.findHomography(pts, pts, error_type="nope")
    with pytest.raises(ValueError, match="Error type"):
        pydegensac.findFundamentalMatrix(pts, pts, error_type="nope")
```

Run: `python -m pytest tests/test_input_validation.py -v`
Expected: `test_failure_path_mask_is_bool_ndarray` FAILS (mask is a list today); the ValueError test may already pass.

- [ ] **Step 2: Fix utils.py**

In `findHomography` (utils.py:106-109): replace `mask = [False]*len(mask)` with `mask = np.zeros(n, dtype=bool)`. In `findFundamentalMatrix` (utils.py:147-149): replace `mask = [False]*n` with `mask = np.zeros(n, dtype=bool)`. Replace the bare `except:` at line 9 with `except ImportError:`, and the two bare `except:` around the error-type dict lookups (lines 93, 133) with `except (KeyError, AttributeError):` keeping the existing `raise ValueError(...)` bodies.

- [ ] **Step 3: Run tests, full gate, commit**

```bash
python -m pytest tests -q
git add src/pydegensac/utils.py tests/test_input_validation.py
git commit -m "fix failure-path mask type and bare excepts in utils.py"
```
Expected: all PASS (goldens exercise the success path — unchanged).

---

### Task 11: Fix bindings bugs (bindings.cpp)

**Files:**
- Modify: `src/pydegensac/bindings.cpp`
- Test: `tests/test_input_validation.py` (extend)

**Interfaces:**
- Consumes: audit findings (uninitialized H/F; missing c_style; no ndim/DIM-match checks; LAF-with-Nx2 OOB; unchecked malloc; deprecated PYBIND11_PLUGIN).
- Produces: bindings that reject malformed input with ValueError and never return uninitialized memory. Success-path outputs bit-identical (golden gate).

- [ ] **Step 1: Write failing tests (append to tests/test_input_validation.py)**

```python
def test_fortran_order_gives_same_result_as_c_order():
    rs = np.random.RandomState(1)
    pts1 = rs.rand(60, 2) * 100
    H_true = np.array([[1.1, 0.05, 3.0], [0.02, 0.95, -2.0], [1e-4, 2e-4, 1.0]])
    p1h = np.concatenate([pts1, np.ones((60, 1))], axis=1)
    p2h = (H_true @ p1h.T).T
    pts2 = p2h[:, :2] / p2h[:, 2:3]
    Hc, mc = pydegensac.findHomography(pts1, pts2, 1.0, seed=7)
    Hf, mf = pydegensac.findHomography(np.asfortranarray(pts1),
                                       np.asfortranarray(pts2), 1.0, seed=7)
    np.testing.assert_array_equal(Hc, Hf)
    np.testing.assert_array_equal(np.asarray(mc), np.asarray(mf))


def test_1d_input_raises():
    with pytest.raises((ValueError, RuntimeError)):
        pydegensac.findHomography_(np.zeros(8), np.zeros(8),
                                   1.0, 0.999, 100, 0, True, 0.0, -1)


def test_mismatched_width_raises():
    with pytest.raises((ValueError, RuntimeError)):
        pydegensac.findHomography_(np.zeros((10, 6)), np.zeros((10, 2)),
                                   1.0, 0.999, 100, 0, True, 0.0, -1)


def test_laf_coef_with_2d_input_raises():
    pts = np.random.RandomState(2).rand(10, 2)
    with pytest.raises((ValueError, RuntimeError)):
        pydegensac.findHomography_(pts, pts, 1.0, 0.999, 100, 0, True, 0.5, -1)
    with pytest.raises((ValueError, RuntimeError)):
        pydegensac.findFundamentalMatrix_(pts, pts, 1.0, 0.999, 100, 0, True,
                                          0.5, True, -1)
```

Run: `python -m pytest tests/test_input_validation.py -v`
Expected: the three `raises` tests FAIL (no exception today — UB instead); the Fortran-order test FAILS (garbage model on F-order input).

- [ ] **Step 2: Fix bindings.cpp**

In BOTH `findHomography_` and `findFundamentalMatrix_`:

```cpp
// 1. signature: force dense C-layout copies of arbitrary input
py::array_t<double, py::array::c_style | py::array::forcecast> x1y1,
py::array_t<double, py::array::c_style | py::array::forcecast> x2y2,

// 2. before reading shape[1]:
if (buf1.ndim != 2 || buf1a.ndim != 2)
    throw std::invalid_argument("points must be 2-D arrays [N,2] or [N,6]");
// 3. after reading DIM/DIMa:
if (DIM != DIMa)
    throw std::invalid_argument("pts1 and pts2 must have the same number of columns");
if (DIM != 2 && DIM != 6)
    throw std::invalid_argument("points must be [N,2] or [N,6]");
// 4. LAF guard:
if (laf_coef > 0 && DIM == 2)
    throw std::invalid_argument("laf_coef > 0 requires [N,6] input");
// 5. zero-init outputs (H at line ~110, F at line ~323):
double H[3*3] = {0};
double F[3*3] = {0};
// 6. every malloc: check and throw std::bad_alloc / std::runtime_error on NULL
```

Also replace `PYBIND11_PLUGIN(pydegensac)` (line ~474) with `PYBIND11_MODULE(pydegensac, m)` (delete the `return m.ptr();` and the module object construction, keep all `m.def(...)` calls unchanged).

- [ ] **Step 3: Rebuild, run all tests, commit**

```bash
pip install . --force-reinstall --no-deps -q && python -m pytest tests -q
git add src/pydegensac/bindings.cpp tests/test_input_validation.py
git commit -m "harden bindings: contiguity, shape checks, zero-init outputs, PYBIND11_MODULE"
```
Expected: everything PASSES including golden bit-equivalence (all fixes are no-ops for valid C-contiguous input; `forcecast` already copied non-double input before, so success paths are byte-identical).

---

### Task 12: Fix build-system bugs (setup.py, CMakeLists.txt)

**Files:**
- Modify: `setup.py`, `CMakeLists.txt`

**Interfaces:**
- Produces: py3.12-safe setup.py; CMake that honors `--debug`; no stale metadata.

- [ ] **Step 1: setup.py**

Delete `from distutils.version import LooseVersion` (line 8) and `copy_test_file` (lines 97-110). Replace the Windows CMake version check (lines 30-35) with:

```python
if platform.system() == "Windows":
    m = re.search(r"version\s*([\d.]+)", out.decode())
    ver = tuple(int(x) for x in m.group(1).split(".")[:3]) if m else (0,)
    if ver < (3, 1):
        raise RuntimeError("CMake >= 3.1.0 is required on Windows")
```

Delete the stale `download_url=` line (127).

- [ ] **Step 2: CMakeLists.txt**

Replace `SET(CMAKE_BUILD_TYPE "RELEASE")` (line 11) with:

```cmake
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
```

(setup.py always passes `-DCMAKE_BUILD_TYPE=Release` for normal builds, so release artifacts are unchanged; `--debug` now actually produces Debug.)

- [ ] **Step 3: Rebuild both configs, gate, commit**

```bash
pip install . --force-reinstall --no-deps -q && python -m pytest tests -q
python setup.py build_ext --debug 2>&1 | grep -i "Build type\|CMAKE_BUILD_TYPE" || true
git add setup.py CMakeLists.txt
git commit -m "fix distutils import, dead code and metadata in setup.py; honor build type in CMake"
```
Expected: release build passes the full suite; the debug build invocation shows Debug config.

---

### Task 13: Run unit tests in CI

**Files:**
- Modify: `.github/workflows/build_wheels.yml` (all three jobs)

**Interfaces:**
- Consumes: `tests/` incl. golden data (committed), `PYDEGENSAC_GOLDEN_EXACT` switch from Task 4.
- Produces: every wheel is tested with the full pytest suite in sanity mode.

- [ ] **Step 1: Extend cibuildwheel test config in each of the three jobs**

For each job's `env:` block change:

```yaml
CIBW_TEST_REQUIRES: opencv-python-headless matplotlib pytest
CIBW_TEST_COMMAND: >
  python -m pytest {project}/tests -q &&
  python -c "import os, runpy; os.chdir(r'{project}/examples');
  runpy.run_path('simple-example.py', run_name='__main__')"
CIBW_ENVIRONMENT: PYDEGENSAC_GOLDEN_EXACT=0
```

(Windows job: merge `PYDEGENSAC_GOLDEN_EXACT=0` into the existing `CIBW_ENVIRONMENT_WINDOWS` line instead of adding a new key. If the machine that captured the goldens matches a CI target toolchain exactly, that job may set `=1` later — start everywhere with sanity mode.)

- [ ] **Step 2: Push a branch and watch one CI run**

```bash
git checkout -b cleanup-golden-tests
git push -u origin cleanup-golden-tests
gh run watch
```
Expected: all three wheel jobs green with pytest output in the logs.

- [ ] **Step 3: Commit is already part of the branch; open the PR**

```bash
gh pr create --title "Golden regression tests, dead-code cleanup, audit bug fixes" --fill
```

---

## Deferred (explicitly NOT in this plan)

- `resids=NULL` switch in bindings (verified pure speed-up: kills per-iteration reallocs+memcpys of an unread buffer) and any live-code dedup/merging — the speed-up sessions. Note the half-guarded exported functions from the audit (`exp_ransacF`, `exp_inFrani`, `exp_inHrani`) are *deleted* by Tasks 6-7, which resolves that finding's latent-crash surface.
- Evaluating the `exp_ranF.h:10` IDEA (adaptive LO triggering) — future work; the comment is preserved for it.
