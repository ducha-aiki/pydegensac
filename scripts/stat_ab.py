#!/usr/bin/env python3
"""Statistical A/B for deliberate numerics changes.

Runs two builds over the golden pairs at many seeds and compares the resulting
distributions. Used to gate changes that alter floating-point results without
altering the algorithm -- the protocol first used for the per-iteration reseed
removal (4430bc7), kept here so it is reproducible rather than ad hoc.

    python3 scripts/stat_ab.py --build-a .ab/pkg-base --build-b .ab/pkg-branch

Each build directory is a `pip install --target` tree containing pydegensac.
The two builds run in separate subprocesses: pydegensac can only be imported
once per process, and the point is to compare two of them.

Metrics per run, following the 2026-08-10 protocol: inlier count, median error
against ground truth over the reported inliers, and the GT precision of those
inliers (the fraction that are GT-consistent at 3 px).
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
DATA_DIR = REPO / "tests" / "data"

#: A p-value below this stops the gate. Not a knob to widen: the 4430bc7 run
#: reported min p = 0.108 over the same 30 comparisons.
FAIL_P = 0.01

# Worker: runs in a subprocess with one build on sys.path, writes JSON to
# stdout. Kept as a string so this stays a single file.
WORKER = r"""
import json, sys
import numpy as np
sys.path.insert(0, sys.argv[1])
import pydegensac

path, runs = sys.argv[2], int(sys.argv[3])
d = np.load(path)
pts1, pts2 = d["pts1"], d["pts2"]
is_h = "H_gt" in d
GT_PX = 3.0


def gt_err(mask):
    # Error of the selected correspondences against ground truth.
    p1 = np.concatenate([pts1[mask], np.ones((int(mask.sum()), 1))], axis=1)
    if is_h:
        proj = (d["H_gt"] @ p1.T).T
        return np.linalg.norm(proj[:, :2] / proj[:, 2:3] - pts2[mask], axis=1)
    p2 = np.concatenate([pts2[mask], np.ones((int(mask.sum()), 1))], axis=1)
    Fx1 = (d["F_gt"] @ p1.T).T
    num = np.abs(np.sum(p2 * Fx1, axis=1))
    return num / np.sqrt(Fx1[:, 0] ** 2 + Fx1[:, 1] ** 2)


gt_inlier = gt_err(np.ones(len(pts1), bool)) <= GT_PX

out = {"inliers": [], "gt_err": [], "gt_prec": []}
for seed in range(1, runs + 1):
    fn = pydegensac.findHomography if is_h else pydegensac.findFundamentalMatrix
    M, mask = fn(pts1, pts2, px_th=float(d["px_th"]), conf=float(d["conf"]),
                 max_iters=int(d["max_iters"]), seed=seed)
    mask = np.asarray(mask, bool)
    n = int(mask.sum())
    out["inliers"].append(n)
    out["gt_err"].append(float(np.median(gt_err(mask))) if n else float("nan"))
    out["gt_prec"].append(float(gt_inlier[mask].mean()) if n else float("nan"))
print(json.dumps(out))
"""


def run_build(build_dir, npz_path, runs):
    proc = subprocess.run(
        [sys.executable, "-c", WORKER, str(build_dir), str(npz_path), str(runs)],
        capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"worker failed for {npz_path.name} with "
                         f"{build_dir}:\n{proc.stderr}")
    return json.loads(proc.stdout)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--build-a", required=True, help="reference build directory")
    ap.add_argument("--build-b", required=True, help="candidate build directory")
    ap.add_argument("--runs", type=int, default=1000, help="seeds per pair")
    args = ap.parse_args()

    pairs = sorted(DATA_DIR.glob("golden_*.npz"))
    if not pairs:
        raise SystemExit(f"no golden fixtures in {DATA_DIR}")

    print(f"{'pair':52s} {'metric':9s} {'KS':>8} {'p':>10}")
    worst = (1.0, None)
    for npz in pairs:
        a = run_build(args.build_a, npz, args.runs)
        b = run_build(args.build_b, npz, args.runs)
        for metric in ("inliers", "gt_err", "gt_prec"):
            xa = np.asarray(a[metric], float)
            xb = np.asarray(b[metric], float)
            xa, xb = xa[np.isfinite(xa)], xb[np.isfinite(xb)]
            ks, p = stats.ks_2samp(xa, xb)
            print(f"{npz.stem[:52]:52s} {metric:9s} {ks:8.4f} {p:10.4g}")
            if p < worst[0]:
                worst = (p, f"{npz.stem} / {metric}")
    print(f"\nmin p = {worst[0]:.4g}  ({worst[1]})   over "
          f"{len(pairs) * 3} comparisons, {args.runs} seeds each")
    if worst[0] < FAIL_P:
        print("FAIL: distributions differ")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
