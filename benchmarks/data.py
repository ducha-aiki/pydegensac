"""Pair iterators for the two benchmark datasets.

Both come from the CVPR-2020 RANSAC tutorial. Fetch them with
``python setup_data.py`` first.

F: one PhotoTourism validation scene. Each pair carries RootSIFT-8k mutual-NN
correspondences, their SNN ratios (``match_conf``, lower = better) and the
ground-truth calibration of both images.

H: EVD and HPatchesSeq. Correspondences come pre-filtered at SNN ~0.85, so
``scores`` are the surviving ratios. Each row of ``matches.h5[key]`` is
``[x1, y1, x2, y2]``; ``Hgt.h5[key]`` maps image-1 pixels to image-2 pixels
(not inverted).
"""
from pathlib import Path

import cv2
import h5py
import numpy as np

DATA = Path(__file__).resolve().parent / "data"

F_SCENES = ("st_peters_square",)
H_DATASETS = ("EVD", "HPatchesSeq")
H_SPLITS = ("val", "test")


def _load_h5(path):
    with h5py.File(path, "r") as f:
        return {k: f[k][()] for k in f.keys()}


def _f_dir(scene):
    d = DATA / "f_data" / scene
    if not d.exists():
        raise FileNotFoundError(f"{d} missing — run `python setup_data.py`")
    return d


def pair_keys_f(scene=F_SCENES[0]):
    """Sorted pair keys, without loading the correspondences (317 MB)."""
    with h5py.File(_f_dir(scene) / "matches.h5", "r") as f:
        return sorted(f.keys())


def iter_pairs_f(scene=F_SCENES[0], keys=None):
    """Yield one dict per pair: name, pts1, pts2 (N,2), scores (N,),
    K1, K2, R1, R2, T1, T2.

    ``keys`` restricts the iteration; correspondences are read per pair rather
    than all at once, so a 600-pair subset costs ~50 MB instead of 400 MB.
    """
    d = _f_dir(scene)
    K1_K2 = _load_h5(d / "K1_K2.h5")
    R = _load_h5(d / "R.h5")
    T = _load_h5(d / "T.h5")
    with h5py.File(d / "matches.h5", "r") as mf, \
            h5py.File(d / "match_conf.h5", "r") as cf:
        for name in (sorted(mf.keys()) if keys is None else keys):
            m = np.asarray(mf[name][()], np.float64)
            id1, id2 = name.split("-")
            yield {
                "name": name,
                "pts1": np.ascontiguousarray(m[:, :2]),
                "pts2": np.ascontiguousarray(m[:, 2:4]),
                "scores": np.asarray(cf[name][()], np.float64).reshape(-1),
                "K1": np.asarray(K1_K2[name][0][0], np.float64),
                "K2": np.asarray(K1_K2[name][0][1], np.float64),
                "R1": np.asarray(R[id1], np.float64),
                "R2": np.asarray(R[id2], np.float64),
                "T1": np.asarray(T[id1], np.float64).reshape(3),
                "T2": np.asarray(T[id2], np.float64).reshape(3),
            }


def _h_image_shapes(dataset, split, key):
    """(h, w) of both images of a pair.

    The metric only needs the shapes, so the images are read header-only where
    OpenCV allows it — but EVD's .png and HPatches' .ppm both need a decode,
    which is why the H setup downloads images at all.
    """
    imgs = DATA / "homography" / dataset / split / "imgs"
    if dataset == "EVD":
        stem = key.split("-")[0]
        p1, p2 = imgs / "1" / f"{stem}.png", imgs / "2" / f"{stem}.png"
    else:
        p1, p2 = imgs / key[:-4] / "1.ppm", imgs / key[:-4] / f"{key[-1]}.ppm"
    out = []
    for p in (p1, p2):
        img = cv2.imread(str(p), cv2.IMREAD_GRAYSCALE)
        if img is None:
            raise FileNotFoundError(f"cannot read {p}")
        out.append(img.shape[:2])
    return out


def iter_pairs_h(dataset, split="test"):
    """Yield one dict per pair: name, pts1, pts2 (N,2), scores (N,),
    H_gt (3,3), shape1, shape2 (h, w)."""
    if dataset not in H_DATASETS:
        raise ValueError(f"unknown dataset {dataset!r}, expected {H_DATASETS}")
    d = DATA / "homography" / dataset / split
    if not (d / "Hgt.h5").exists():
        raise FileNotFoundError(
            f"{d / 'Hgt.h5'} missing — run `python setup_data.py` "
            f"(the test split's ground truth is recovered from upstream)")
    matches = _load_h5(d / "matches.h5")
    conf = _load_h5(d / "match_conf.h5")
    hgt = _load_h5(d / "Hgt.h5")
    for key in sorted(hgt):
        m = np.asarray(matches[key], np.float64)
        shape1, shape2 = _h_image_shapes(dataset, split, key)
        yield {
            "name": key,
            "pts1": np.ascontiguousarray(m[:, :2]),
            "pts2": np.ascontiguousarray(m[:, 2:4]),
            "scores": np.asarray(conf[key], np.float64).reshape(-1),
            "H_gt": np.asarray(hgt[key], np.float64),
            "shape1": shape1,
            "shape2": shape2,
        }


def split_keys(keys, n_tune, n_eval, seed):
    """Disjoint (tune, eval) key subsets drawn from one seeded permutation."""
    if n_tune + n_eval > len(keys):
        raise ValueError(f"asked for {n_tune}+{n_eval} pairs, have {len(keys)}")
    idx = np.random.default_rng(seed).permutation(len(keys))
    take = lambda sl: [keys[i] for i in sorted(idx[sl])]  # noqa: E731
    return take(slice(0, n_tune)), take(slice(n_tune, n_tune + n_eval))
