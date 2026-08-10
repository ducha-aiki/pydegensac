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


def _h_gt_error(d, mask):
    p1 = np.concatenate([d["pts1"][mask], np.ones((mask.sum(), 1))], axis=1)
    proj = (d["H_gt"] @ p1.T).T
    err = np.linalg.norm(proj[:, :2] / proj[:, 2:3] - d["pts2"][mask], axis=1)
    return np.median(err)


def _f_gt_error(d, mask):
    p1 = np.concatenate([d["pts1"][mask], np.ones((mask.sum(), 1))], axis=1)
    p2 = np.concatenate([d["pts2"][mask], np.ones((mask.sum(), 1))], axis=1)
    Fx1 = (d["F_gt"] @ p1.T).T
    num = np.abs(np.sum(p2 * Fx1, axis=1))
    den = np.sqrt(Fx1[:, 0] ** 2 + Fx1[:, 1] ** 2)
    return np.median(num / den)


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
        # GT-agreement sanity: the freshly computed model+mask must still be
        # geometrically consistent with ground truth, not just internally
        # reproducible.
        if mask.sum() >= 10:
            if "H_gt" in d:
                assert _h_gt_error(d, mask) < 10.0
            else:
                assert _f_gt_error(d, mask) < 5.0


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
