"""Estimator-only timing driver: pydegensac and nothing else.

    PYTHONPATH=<pkgdir> taskset -c 2 python iso_bench.py h|f [--reps N]

Two traps this exists to avoid. Profiling ``run.py`` shows numpy + cv2 at 85%
of process time — that is the reprojection metric, not the estimator. And
unpinned OpenBLAS burns 20-43% in ``__sched_yield`` spinning on 9x9 matrices,
which swamps everything being measured; hence the thread caps below and
``taskset`` in the caller.

Runs with a fixed seed so every build does *identical* work: the printed
inlier checksum must match across arms, and if it does, the only thing that can
differ is how fast that work executes.
"""
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import argparse  # noqa: E402
import json  # noqa: E402
import statistics  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

import pydegensac  # noqa: E402

#: pydegensac's own rows from tuned_config.json, at the top of each budget
#: ladder in methods.py. F is capped at 10000 rather than 50000 to keep a
#: sweep to a few minutes; the ranking does not move with the budget.
CFG = {
    "h": dict(px_th=4.0, conf=0.999, max_iters=25000),
    "f": dict(px_th=0.5, conf=0.9999, max_iters=10000),
}
SEED = 12345
DEFAULT_NPZ = Path(__file__).resolve().parents[1] / ".ab" / "iso" / "pairs.npz"


def load(problem, npz):
    pairs, i = [], 0
    while f"{problem}/{i}/1" in npz:
        pairs.append((npz[f"{problem}/{i}/1"], npz[f"{problem}/{i}/2"]))
        i += 1
    if not pairs:
        raise SystemExit(f"no {problem} pairs in the npz")
    return pairs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("problem", choices=("h", "f"))
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--label", default="")
    ap.add_argument("--npz", default=str(DEFAULT_NPZ))
    a = ap.parse_args()

    if not Path(a.npz).exists():
        raise SystemExit(f"{a.npz} missing — run `python dump_pairs.py` first")

    cfg = CFG[a.problem]
    pairs = load(a.problem, np.load(a.npz))
    est = (pydegensac.findHomography if a.problem == "h"
           else pydegensac.findFundamentalMatrix)

    def one_pass():
        n_in = 0
        t0 = time.perf_counter()
        for p1, p2 in pairs:
            M, mask = est(p1, p2, cfg["px_th"], conf=cfg["conf"],
                          max_iters=cfg["max_iters"], seed=SEED)
            n_in += int(mask.sum())
        return time.perf_counter() - t0, n_in

    one_pass()  # warm up: the first call faults in the pages
    times, checks = [], set()
    for _ in range(a.reps):
        dt, n_in = one_pass()
        times.append(dt)
        checks.add(n_in)
        print(json.dumps({"label": a.label, "s": round(dt, 4),
                          "ms_per_pair": round(1e3 * dt / len(pairs), 3),
                          "inliers": n_in}), flush=True)

    print(json.dumps({
        "label": a.label, "problem": a.problem, "pairs": len(pairs),
        "median_s": round(statistics.median(times), 4),
        "min_s": round(min(times), 4),
        "ms_per_pair": round(1e3 * statistics.median(times) / len(pairs), 3),
        "inlier_checksum": sorted(checks),
        "so": pydegensac.__file__,
    }), flush=True)
    if len(checks) > 1:
        print("WARNING: nondeterministic inlier count across repeats",
              file=sys.stderr)


main()
