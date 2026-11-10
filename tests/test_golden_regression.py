import hashlib
import os
from pathlib import Path

import numpy as np
import pytest

import pydegensac

DATA_DIR = Path(__file__).parent / "data"
H_FILES = sorted(DATA_DIR.glob("golden_h_*.npz"))
F_FILES = sorted(DATA_DIR.glob("golden_f_*.npz"))
SEEDS = (42, 2024)

# Exact fixtures below cover the fixed-size eigensolver and rank-2 projection on the capturing
# macOS/arm64 toolchain. CI uses PYDEGENSAC_GOLDEN_EXACT=0 for other toolchains;
# the NPZ files retain the input correspondences and geometric ground truth.
EXACT = os.environ.get("PYDEGENSAC_GOLDEN_EXACT", "1") == "1"

NUMERICAL_KERNEL_DIGESTS = {('golden_f_reichstag_05461164_9050854768_05466646_5360480312', 'F_seed2024'): '941f4374bc4b34934a2b9f8a74b7f7ca9d7dc771175cbc80c6221ebfd04c6ed4',
 ('golden_f_reichstag_05461164_9050854768_05466646_5360480312', 'F_seed42'): '038958fc7c36602cfa511a7f225f67312d4e81709dbd1aa827f3fb88c595c1ee',
 ('golden_f_reichstag_05534141_6340060522_05545431_1341290848', 'F_seed2024'): 'f7ca926ad8dca6684cc2f778a1494b6ba82890ad68f25a8d0eb3c7e240e13b73',
 ('golden_f_reichstag_05534141_6340060522_05545431_1341290848', 'F_seed42'): '0bb9c462cba75526de2e30c42da0ee49df55a466aecac853e1666091ed8f0f92',
 ('golden_f_reichstag_05545431_1341290848_05791347_12791964625', 'F_seed2024'): '608388d709d8a026e725a73114bd36861af63295b4ed7510c625c723f9f097b5',
 ('golden_f_reichstag_05545431_1341290848_05791347_12791964625', 'F_seed42'): '507f46599d30b63c242ebd8abb86831e0316a8213620d6c792f3b183cf2def9f',
 ('golden_f_reichstag_06229406_8584869180_06373813_9127207861', 'F_seed2024'): 'b70f45be1b04a16220e6e747393d31c96e2b5fed0dd3aa0df7102755dcdd37e6',
 ('golden_f_reichstag_06229406_8584869180_06373813_9127207861', 'F_seed42'): '848d55d53e2f6fc5ba8c443f55e8f931b31150a125641865058a19e1617a317d',
 ('golden_f_reichstag_06373813_9127207861_06639257_4979506057', 'F_seed2024'): '1fa90e33cd8ed18b9223d3ac68288272381d66c639740c915924333ae70d85bd',
 ('golden_f_reichstag_06373813_9127207861_06639257_4979506057', 'F_seed42'): '4fca0f43c84e45e24a36cd3d2475a0cc2e559e78f84b0c6156bb9dadc4d5b17d',
 ('golden_h_adam', 'H_seed2024'): 'd8e51fb4a21545091dead72aecfd78f256b665e86d586ec3862719a57c219d49',
 ('golden_h_adam', 'H_seed42'): 'a1c7124f9f30aae62893530ca7b5e0ceb09af3400562902c5ab3a84ece04ec19',
 ('golden_h_cat', 'H_seed2024'): 'a1e7b9a4ddcfd8ccdcf464ca4afc442dafdc6218b5180920342725a0e9c66ebe',
 ('golden_h_cat', 'H_seed42'): 'ce9afff065ba47fcfd87503af56e1839d1f3c338a1ff1fcb59d959eb6e146114',
 ('golden_h_dum', 'H_seed2024'): '523079de8f144cf77d32f6623ee31fd910493a243c98e549e335be7920fa193c',
 ('golden_h_dum', 'H_seed42'): 'bbb066e90d00594fd6d17b169e3d4295a39a84eeac4fa235efe39b5e76cad927',
 ('golden_h_face', 'H_seed2024'): 'd12ae8cef88a622f750048807091c33684ff8c1dd4a2829efa142fd6fe9b2596',
 ('golden_h_face', 'H_seed42'): '3f1aac4f0a8951af2759557b50f8b41d10f0ca8e63878f420d08ee647194ab52',
 ('golden_h_fox', 'H_seed2024'): '622105911f7a35129a4bc7357a07c424b59d921e662e921ea0672b06f096cfa6',
 ('golden_h_fox', 'H_seed42'): '622105911f7a35129a4bc7357a07c424b59d921e662e921ea0672b06f096cfa6'}

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


def _assert_matches_golden(model, mask, d, model_key, mask_key, rerun, fixture):
    mask = np.asarray(mask, dtype=bool)
    if EXACT:
        digest = hashlib.sha256(model.tobytes() + mask.tobytes()).hexdigest()
        assert digest == NUMERICAL_KERNEL_DIGESTS[(fixture, model_key)]
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
    _assert_matches_golden(H, mask, d, f"H_seed{seed}", f"mask_seed{seed}", call, path.stem)


@pytest.mark.parametrize("path", F_FILES, ids=lambda p: p.stem)
@pytest.mark.parametrize("seed", SEEDS)
def test_fundamental_bit_equivalent(path, seed):
    d = np.load(path)
    call = lambda: pydegensac.findFundamentalMatrix(
        d["pts1"], d["pts2"],
        px_th=float(d["px_th"]), conf=float(d["conf"]),
        max_iters=int(d["max_iters"]), seed=seed)
    F, mask = call()
    _assert_matches_golden(F, mask, d, f"F_seed{seed}", f"mask_seed{seed}", call, path.stem)


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
