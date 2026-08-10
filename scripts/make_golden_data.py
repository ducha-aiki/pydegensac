#!/usr/bin/env python3
"""Generate golden regression data for pydegensac.

The H-side goldens (EVD) use pre-extracted MODS tentative correspondences
that ship with the dataset -- plain SIFT+ratio matching does not find
enough GT-consistent correspondences on EVD's extreme-viewpoint pairs (see
task-2-report.md for the diagnostic). The F-side goldens (Task 3, IMC
pairs) still use imc-2021-simple's SIFT extractor + ratio matching, which
is why those helpers stay in this module. Ground truth comes with the
datasets (EVD homographies, IMC calibration); the "golden" outputs are the
current pydegensac results at fixed seeds. Regenerate ONLY deliberately:
committing new goldens redefines the bit-equivalence baseline.
"""
import os

# sift_eval (below) pulls in imc2021, which imports pycolmap; pycolmap
# bundles its own OpenMP runtime that collides with the one already loaded
# by torch/opencv on macOS ("OMP: Error #15" -> SIGABRT on import). This
# must be set before cv2/torch/pycolmap are imported anywhere in-process.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import subprocess
import sys
import tempfile
import urllib.request
import zipfile
from pathlib import Path

import cv2
import numpy as np

REPO = Path(__file__).resolve().parent.parent
DATA_DIR = REPO / "tests" / "data"
CACHE_DIR = REPO / ".cache"
IMC_SIMPLE = Path("/Users/oldufo/dev/imc-2021-simple")
sys.path.insert(0, str(IMC_SIMPLE / "examples"))
sys.path.insert(0, str(IMC_SIMPLE))
from sift_eval import extract_sift, _find_image as find_image  # noqa: E402

import pydegensac  # noqa: E402

N_FEATURES = 8000
RATIO_TH = 0.8
SEEDS = (42, 2024)
H_PARAMS = dict(px_th=3.0, conf=0.999, max_iters=20000)
F_PARAMS = dict(px_th=0.75, conf=0.9999, max_iters=10000)

EVD_TENTATIVES_URL = "http://cmp.felk.cvut.cz/wbs/datasets/EVD_tentatives.zip"

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


def ensure_evd_tentatives(cache_dir=CACHE_DIR):
    """Download (if needed) and unpack the pre-extracted MODS tentative
    correspondences for EVD. Format, per inspection of one .txt file:
    comma-separated, header row
    ``x1,y1,x2,y2,FGINN_ratio,SNN_ratio,detector,descriptor,is_correct``;
    we only use the x1,y1,x2,y2 columns (image-1 and image-2 pixel coords)
    and deliberately ignore `is_correct` -- pydegensac's RANSAC, not the
    MODS annotation, is what decides inliers here."""
    tent_dir = cache_dir / "EVD_tentatives"
    if not tent_dir.exists():
        cache_dir.mkdir(parents=True, exist_ok=True)
        zip_path = cache_dir / "EVD_tentatives.zip"
        if not zip_path.exists():
            urllib.request.urlretrieve(EVD_TENTATIVES_URL, zip_path)
        with zipfile.ZipFile(zip_path) as z:
            z.extractall(tent_dir)
    return tent_dir


def load_evd_tentatives(tent_dir, name):
    """Load MODS tentative correspondences for one EVD pair by name."""
    fname = tent_dir / f"{name}.png_m.txt"
    if not fname.exists():
        return None, None
    data = np.genfromtxt(fname, delimiter=",", skip_header=1, usecols=(0, 1, 2, 3))
    if data.ndim == 1:
        data = data.reshape(1, -1)
    pts1 = data[:, 0:2].astype(np.float64)
    pts2 = data[:, 2:4].astype(np.float64)
    return pts1, pts2


def make_evd_goldens(max_pairs=5):
    from wxbs_benchmark.dataset import EVDDataset
    dset = EVDDataset(str(REPO / ".cache" / "EVD"), download=True)
    tent_dir = ensure_evd_tentatives()
    saved = 0
    for pair in dset:
        keys = set(pair.keys())
        h_key = "H" if "H" in keys else ("H_gt" if "H_gt" in keys else None)
        assert h_key is not None, f"no homography key in EVD sample: {keys}"
        name = pair.get("name", f"pair{saved}")
        pts1, pts2 = load_evd_tentatives(tent_dir, name)
        if pts1 is None or len(pts1) < 20:
            print(f"skip {name}: too few tentative correspondences")
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
    assert saved >= 3, f"only {saved} EVD pairs qualified out of {len(dset)}"


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
    # determinism check: capture must be reproducible before it is golden
    F2, mask2 = pydegensac.findFundamentalMatrix(pts1, pts2, seed=seed, **F_PARAMS)
    assert np.array_equal(F, F2) and np.array_equal(mask, np.asarray(mask2, dtype=bool))
    return F, mask


def f_sanity(F, mask, pts1, pts2, c1, c2, max_deg=15.0):
    """Pose recovered from F must agree with GT calibration within max_deg."""
    from imc2021.metrics import pose_error
    if mask.sum() < 30:
        return False
    err = pose_error(F, pts1[mask], pts2[mask], c1, c2)
    return err["max_err"] < max_deg


def _ensure_phototourism_extracted(dest, dataset="phototourism", split="val"):
    """Work around a py311 incompatibility in imc2021.download.download_dataset:
    it calls ``TarFile.extractall(path=dest, filter="fully_trusted")``, but the
    ``filter`` kwarg was only added to extractall in Python 3.12 (PEP 706) --
    this raises TypeError on the active py311 env. If the tarball is already
    downloaded (as ensured by this script) but not yet extracted, extract it
    ourselves and write the same sentinel file download_dataset() checks
    (same path, same "<scene>\\n..." content), so its own call becomes a
    no-op (sentinel present -> returns immediately, buggy extractall never
    runs). No changes to the imc-2021-simple package itself."""
    import tarfile
    from imc2021.download import _CHECKSUMS, _scene_names_from_tar

    dest = Path(dest)
    out_dir = dest / dataset
    sentinel = out_dir / f".{split}.done"
    if sentinel.exists():
        return
    filename, _ = _CHECKSUMS[(dataset, split)]
    tarball_path = dest / filename
    if not tarball_path.exists():
        return  # not downloaded yet -- let download_dataset() do its normal thing
    print(f"Extracting {tarball_path} (py311 workaround for extractall(filter=...))")
    with tarfile.open(tarball_path, "r:gz") as tar:
        scene_names = _scene_names_from_tar(tar, dataset)
        tar.extractall(path=dest)
    if not out_dir.exists():
        raise RuntimeError(f"Extraction did not produce expected directory {out_dir}")
    sentinel.write_text("\n".join(scene_names) + "\n")


def make_imc_goldens(max_pairs=5, scene="reichstag"):
    from imc2021.download import download_dataset
    from imc2021.io import load_calibration
    imc_dest = REPO / ".cache" / "imc"
    _ensure_phototourism_extracted(imc_dest)
    root = Path(download_dataset("phototourism", "val", dest=str(imc_dest)))
    scene_dir = root / scene / "set_100"
    calib = load_calibration(scene_dir)
    names = sorted(calib.keys())
    images_dir = scene_dir / "images"
    feats = {}
    saved = 0
    # walk consecutive-ish pairs until enough qualify
    for i in range(0, len(names) - 1):
        a, b = names[i], names[i + 1]
        for n in (a, b):
            if n not in feats:
                img_path = find_image(images_dir, n)
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
            print(f"skip {a}-{b}: sanity check failed")
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


if __name__ == "__main__":
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    make_evd_goldens()
    make_imc_goldens()
