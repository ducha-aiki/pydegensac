"""Run the benchmark. Three subcommands: ``tune``, ``f``, ``h``.

    python run.py tune f            # pick thresholds on the F tuning subset
    python run.py tune h            # ... on the H val split
    python run.py f                 # F time-mAA sweep -> results/f.jsonl
    python run.py h                 # H time-mAA sweep -> results/h.jsonl

Timing is single-threaded (``OMP_NUM_THREADS=1`` is set before numpy loads,
``cv2.setNumThreads(1)`` after) and covers the estimator call only. Methods are
interleaved per pair so machine drift hits all of them equally, and each method
is warmed up before the timed loop.
"""
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import argparse  # noqa: E402
import itertools  # noqa: E402
import json  # noqa: E402
import platform  # noqa: E402
import subprocess  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import cv2  # noqa: E402
import numpy as np  # noqa: E402

import data  # noqa: E402
import methods  # noqa: E402
import metrics  # noqa: E402

cv2.setNumThreads(1)

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
CONFIG_PATH = HERE / "tuned_config.json"

#: Tuning grids, matching the reference studies.
PX_GRID = (0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 4.0)
RATIO_GRID = (0.6, 0.65, 0.7, 0.75, 0.8, 0.85)

#: Pair counts for the F scene: disjoint tuning and evaluation subsets.
F_N_TUNE = 300
F_N_EVAL = 600
F_SEED = 0


def load_config():
    return json.loads(CONFIG_PATH.read_text())


def git_sha():
    try:
        return subprocess.check_output(
            ["git", "-C", str(HERE), "rev-parse", "--short", "HEAD"],
            text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return "unknown"


def meta(label, extra):
    return {"record": "meta", "label": label, "git_sha": git_sha(),
            "host": platform.node(), "python": platform.python_version(),
            "versions": methods.versions(), **extra}


# --------------------------------------------------------------------------
# Per-pair evaluation
# --------------------------------------------------------------------------

def _filter(pair, ratio_th):
    """Apply the SNN ratio test. ``None`` keeps everything (the H data is
    already filtered upstream and has no ratio knob)."""
    if ratio_th is None or ratio_th >= 1.0:
        return pair["pts1"], pair["pts2"], pair["scores"]
    keep = pair["scores"] <= ratio_th
    return pair["pts1"][keep], pair["pts2"][keep], pair["scores"][keep]


def _score_f(pair, F, mask, pts1, pts2):
    if F is None or mask.sum() < 5:
        return metrics.FAIL_ERR
    return metrics.pose_error(F, pts1[mask], pts2[mask], pair["K1"], pair["K2"],
                              pair["R1"], pair["R2"], pair["T1"], pair["T2"])


def _score_h(pair, H, mask, pts1, pts2):
    if H is None or mask.sum() < 4:
        return metrics.FAIL_ERR
    return metrics.reprojection_error(pair["shape1"], pair["shape2"],
                                      pair["H_gt"], H)


def evaluate(problem, pair, fn, px_th, ratio_th, max_iters):
    """One timed estimator call -> (error, seconds, n_inliers)."""
    pts1, pts2, scores = _filter(pair, ratio_th)
    t0 = time.perf_counter()
    try:
        M, mask = fn(pts1, pts2, px_th, scores, max_iters=max_iters)
    except Exception as exc:  # a backend blowing up is a failure, not a crash
        print(f"  ! {pair['name']}: {type(exc).__name__}: {exc}")
        M, mask = None, np.zeros(len(pts1), bool)
    dt = time.perf_counter() - t0
    score = _score_f if problem == "f" else _score_h
    return score(pair, M, mask, pts1, pts2), dt, int(mask.sum())


def warmup(problem, fn, pair, px_th, ratio_th):
    """One untimed call, so first-call costs don't land on the first pair."""
    try:
        evaluate(problem, pair, fn, px_th, ratio_th, 100)
    except Exception:
        pass


# --------------------------------------------------------------------------
# tune
# --------------------------------------------------------------------------

def cmd_tune(args):
    problem = args.problem
    registry = methods.registry(problem)
    names = args.methods or list(methods.METHOD_NAMES)
    top_budget = methods.budgets(problem)[-1]

    if problem == "f":
        keys, _ = data.split_keys(data.pair_keys_f(), F_N_TUNE, F_N_EVAL,
                                  F_SEED)
        pairs = list(data.iter_pairs_f(keys=keys))
        grid = list(itertools.product(PX_GRID, RATIO_GRID))
        subsets = {"st_peters_square": pairs}
    else:
        grid = [(px, None) for px in PX_GRID]
        subsets = {ds: list(data.iter_pairs_h(ds, "val"))
                   for ds in data.H_DATASETS}

    print(f"tuning {problem} on "
          + ", ".join(f"{k}: {len(v)} pairs" for k, v in subsets.items())
          + f" at max_iters={top_budget}")

    out = {}
    for name in names:
        fn = registry[name]
        best = None
        for px_th, ratio_th in grid:
            scores = {}
            for subset, pairs in subsets.items():
                warmup(problem, fn, pairs[0], px_th, ratio_th)
                errs = [evaluate(problem, p, fn, px_th, ratio_th, top_budget)[0]
                        for p in pairs]
                scores[subset] = (metrics.maa_f(errs) if problem == "f"
                                  else metrics.maa_h(errs))
            # One config per method: rank on the mean over subsets so a method
            # cannot be tuned separately per dataset while others are not.
            mean = float(np.mean(list(scores.values())))
            print(f"  {name:16s} px={px_th:<5} ratio={ratio_th} "
                  + " ".join(f"{k}={v:.4f}" for k, v in scores.items())
                  + f" mean={mean:.4f}", flush=True)
            if best is None or mean > best["maa"]:
                best = {"px_th": px_th, "maa": mean, "per_subset": scores}
                if ratio_th is not None:
                    best["ratio_th"] = ratio_th
        out[name] = best
        print(f"  -> {name}: {best}", flush=True)

    dest = RESULTS / f"tuning_{problem}.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(
        {"meta": meta(args.label, {"problem": problem, "grid": grid,
                                   "max_iters": top_budget}),
         "best": out}, indent=1))
    print(f"wrote {dest}")


# --------------------------------------------------------------------------
# sweeps
# --------------------------------------------------------------------------

def _sweep(problem, subsets, names, config, budgets, out_path, label):
    registry = methods.registry(problem)
    missing = [n for n in names if config.get(n, {}).get("px_th") is None]
    if missing:
        raise SystemExit(
            f"no tuned threshold for {missing} in {CONFIG_PATH.name} — "
            f"run `python run.py tune {problem}` first")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_total = sum(len(p) for p in subsets.values())
    with open(out_path, "w") as fh:
        fh.write(json.dumps(meta(label, {
            "problem": problem, "budgets": list(budgets),
            "config": {n: config[n] for n in names},
            "subsets": {k: len(v) for k, v in subsets.items()}})) + "\n")

        plan = [(name, registry[name], config[name]["px_th"],
                 config[name].get("ratio_th")) for name in names]
        for name, fn, px_th, ratio_th in plan:
            for pairs in subsets.values():
                warmup(problem, fn, pairs[0], px_th, ratio_th)

        done = 0
        t_start = time.perf_counter()
        # Pair-major: every method sees the same pair back to back, so machine
        # drift over the run cannot favour whichever method ran first.
        for subset, pairs in subsets.items():
            for pair in pairs:
                for name, fn, px_th, ratio_th in plan:
                    for budget in budgets:
                        err, dt, n_in = evaluate(problem, pair, fn, px_th,
                                                 ratio_th, budget)
                        fh.write(json.dumps({
                            "record": "pair", "method": name, "subset": subset,
                            "pair": pair["name"], "budget": budget,
                            "px_th": px_th, "ratio_th": ratio_th,
                            "err": None if not np.isfinite(err) else err,
                            "time": dt, "n_inliers": n_in}) + "\n")
                done += 1
                if done % 25 == 0 or done == n_total:
                    el = time.perf_counter() - t_start
                    print(f"  {done}/{n_total} pairs ({el:.0f}s)", flush=True)
                    fh.flush()
    print(f"wrote {out_path}")


def cmd_f(args):
    config = load_config()["f"]
    names = args.methods or list(methods.METHOD_NAMES)
    _, keys = data.split_keys(data.pair_keys_f(), F_N_TUNE, F_N_EVAL, F_SEED)
    if args.limit:
        keys = keys[:args.limit]
    pairs = list(data.iter_pairs_f(keys=keys))
    print(f"F sweep: {len(pairs)} pairs, methods {names}")
    _sweep("f", {"st_peters_square": pairs}, names, config,
           methods.F_BUDGETS, RESULTS / args.out, args.label)


def cmd_h(args):
    config = load_config()["h"]
    names = args.methods or list(methods.METHOD_NAMES)
    subsets = {}
    for ds in data.H_DATASETS:
        pairs = list(data.iter_pairs_h(ds, args.split))
        subsets[ds] = pairs[:args.limit] if args.limit else pairs
    print(f"H sweep ({args.split}): "
          + ", ".join(f"{k}: {len(v)} pairs" for k, v in subsets.items())
          + f", methods {names}")
    _sweep("h", subsets, names, config, methods.H_BUDGETS,
           RESULTS / args.out, args.label)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p, default_out):
        p.add_argument("--methods", nargs="*", choices=methods.METHOD_NAMES,
                       help="default: all")
        p.add_argument("--label", default="local",
                       help="tag for this run, e.g. the pydegensac build")
        p.add_argument("--out", default=default_out)
        p.add_argument("--limit", type=int, default=0,
                       help="cap pairs per subset (smoke tests)")

    t = sub.add_parser("tune", help="grid-search thresholds")
    t.add_argument("problem", choices=["f", "h"])
    t.add_argument("--methods", nargs="*", choices=methods.METHOD_NAMES)
    t.add_argument("--label", default="local")
    t.set_defaults(func=cmd_tune)

    f = sub.add_parser("f", help="fundamental-matrix time-mAA sweep")
    common(f, "f.jsonl")
    f.set_defaults(func=cmd_f)

    h = sub.add_parser("h", help="homography time-mAA sweep")
    common(h, "h.jsonl")
    h.add_argument("--split", default="test", choices=data.H_SPLITS)
    h.set_defaults(func=cmd_h)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
