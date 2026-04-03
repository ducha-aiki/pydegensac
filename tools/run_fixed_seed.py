#!/usr/bin/env python
import argparse
import json
from pathlib import Path
from time import perf_counter

import cv2
import numpy as np

import pydegensac


ROOT = Path(__file__).resolve().parents[1]
IMG_DIR = ROOT / "examples" / "img" / "v_dogman"

H_TH = 4.0
H_CONF = 0.99
H_ITERS = 2000

F_TH = 0.5
F_CONF = 0.999
F_ITERS = 50000


def load_tentatives():
    img1 = cv2.cvtColor(cv2.imread(str(IMG_DIR / "1.ppm")), cv2.COLOR_BGR2RGB)
    img2 = cv2.cvtColor(cv2.imread(str(IMG_DIR / "6.ppm")), cv2.COLOR_BGR2RGB)

    det = cv2.AKAZE_create(descriptor_type=3, threshold=0.00001)
    kps1, descs1 = det.detectAndCompute(img1, None)
    kps2, descs2 = det.detectAndCompute(img2, None)

    bf = cv2.BFMatcher()
    matches = bf.knnMatch(descs1, descs2, k=2)

    tentatives = []
    for first, second in matches:
        if first.distance < 0.9 * second.distance:
            tentatives.append(first)

    src_pts = np.float64([kps1[m.queryIdx].pt for m in tentatives]).reshape(-1, 2)
    dst_pts = np.float64([kps2[m.trainIdx].pt for m in tentatives]).reshape(-1, 2)
    return src_pts, dst_pts


def matrix_signature(matrix):
    return " ".join(f"{value:.17g}" for value in matrix.reshape(-1))


def run_estimation(seed):
    src_pts, dst_pts = load_tentatives()

    start = perf_counter()
    H, H_mask = pydegensac.findHomography(
        src_pts,
        dst_pts,
        H_TH,
        H_CONF,
        H_ITERS,
        seed=seed,
    )
    h_time = perf_counter() - start

    start = perf_counter()
    F, F_mask = pydegensac.findFundamentalMatrix(
        src_pts,
        dst_pts,
        F_TH,
        F_CONF,
        F_ITERS,
        enable_degeneracy_check=True,
        seed=seed,
    )
    f_time = perf_counter() - start

    H_mask = np.asarray(H_mask, dtype=bool)
    F_mask = np.asarray(F_mask, dtype=bool)

    return {
        "seed": seed,
        "num_tentatives": int(src_pts.shape[0]),
        "homography": {
            "runtime_sec": h_time,
            "inliers": int(H_mask.sum()),
            "matrix": np.asarray(H, dtype=np.float64).tolist(),
            "mask": H_mask.astype(int).tolist(),
            "signature": matrix_signature(np.asarray(H, dtype=np.float64)),
        },
        "fundamental": {
            "runtime_sec": f_time,
            "inliers": int(F_mask.sum()),
            "matrix": np.asarray(F, dtype=np.float64).tolist(),
            "mask": F_mask.astype(int).tolist(),
            "signature": matrix_signature(np.asarray(F, dtype=np.float64)),
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    result = run_estimation(args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
