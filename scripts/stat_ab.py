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

Metrics per run. The first three follow the 2026-08-10 protocol: inlier count,
median error against ground truth over the reported inliers, and the GT
precision of those inliers (the fraction GT-consistent at 3 px). All three are
functions of the inlier *mask* alone, so they are blind to a change that
perturbs the returned model without moving the inlier set -- which is the usual
effect of a rounding change. `model_err` closes that hole by scoring the
returned model itself against ground truth.
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
#: reported min p = 0.108 over the mask-based metrics.
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


def model_err(M):
    # Quality of the returned MODEL against ground truth, measured on the
    # GT-consistent correspondences. The three mask-based metrics above are
    # blind to a model whose last bits move without the inlier set changing,
    # which is exactly what a rounding change usually does.
    if M is None or not np.any(M):
        return float("nan")
    sel = gt_inlier
    if sel.sum() < 5:
        return float("nan")
    p1 = np.concatenate([pts1[sel], np.ones((int(sel.sum()), 1))], axis=1)
    if is_h:
        proj = (M @ p1.T).T
        w = proj[:, 2:3]
        w = np.where(np.abs(w) < 1e-12, 1e-12, w)
        ref = (d["H_gt"] @ p1.T).T
        wr = ref[:, 2:3]
        wr = np.where(np.abs(wr) < 1e-12, 1e-12, wr)
        return float(np.median(np.linalg.norm(
            proj[:, :2] / w - ref[:, :2] / wr, axis=1)))
    # F: Sampson error of the estimated model on GT-consistent matches
    p2 = np.concatenate([pts2[sel], np.ones((int(sel.sum()), 1))], axis=1)
    Fx1 = (M @ p1.T).T
    Ftx2 = (M.T @ p2.T).T
    num = np.sum(p2 * Fx1, axis=1) ** 2
    den = Fx1[:, 0]**2 + Fx1[:, 1]**2 + Ftx2[:, 0]**2 + Ftx2[:, 1]**2
    den = np.where(den < 1e-12, 1e-12, den)
    return float(np.median(num / den))


out = {"inliers": [], "gt_err": [], "gt_prec": [], "model_err": []}
for seed in range(1, runs + 1):
    fn = pydegensac.findHomography if is_h else pydegensac.findFundamentalMatrix
    M, mask = fn(pts1, pts2, px_th=float(d["px_th"]), conf=float(d["conf"]),
                 max_iters=int(d["max_iters"]), seed=seed)
    mask = np.asarray(mask, bool)
    n = int(mask.sum())
    out["inliers"].append(n)
    out["gt_err"].append(float(np.median(gt_err(mask))) if n else float("nan"))
    out["gt_prec"].append(float(gt_inlier[mask].mean()) if n else float("nan"))
    out["model_err"].append(model_err(np.asarray(M, float) if M is not None else None))
print(json.dumps(out))
"""


#: Significant digits kept before the KS test.
#:
#: These estimators return a handful of distinct outcomes across seeds, so the
#: samples are heavily tied. A rounding change moves every tied value by ~1e-13
#: without changing which outcome was found, and a two-sample KS test reads
#: that as the CDFs stepping at different places -- the statistic jumps to the
#: height of the tie (measured: KS 0.86, p = 0, on distributions whose medians
#: and means agreed to four decimals). Quantising first makes the test see
#: "same outcome" rather than "different by one ulp". Nine digits is far below
#: any difference that would matter: model errors here are 0.05-2 px.
QUANT_DIGITS = 9


def _quantize(x):
    if x.size == 0:
        return x
    scale = np.maximum(np.abs(x), 1e-300)
    exp = np.floor(np.log10(scale))
    step = 10.0 ** (exp - (QUANT_DIGITS - 1))
    return np.round(x / step) * step


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
        for metric in ("inliers", "gt_err", "gt_prec", "model_err"):
            xa = np.asarray(a[metric], float)
            xb = np.asarray(b[metric], float)
            xa, xb = xa[np.isfinite(xa)], xb[np.isfinite(xb)]
            xa, xb = _quantize(xa), _quantize(xb)
            ks, p = stats.ks_2samp(xa, xb)
            print(f"{npz.stem[:52]:52s} {metric:9s} {ks:8.4f} {p:10.4g}")
            if p < worst[0]:
                worst = (p, f"{npz.stem} / {metric}")
    print(f"\nmin p = {worst[0]:.4g}  ({worst[1]})   over "
          f"{len(pairs) * 4} comparisons, {args.runs} seeds each")
    if worst[0] < FAIL_P:
        print("FAIL: distributions differ")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
