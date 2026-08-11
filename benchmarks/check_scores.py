"""Check which way each dataset's per-match score points, against ground truth.

PROSAC needs correspondences ordered best-first, so the benchmark has to know
whether a dataset's ``match_conf`` is "lower is better" (an SNN ratio) or the
opposite. The CVPR-2020 RANSAC tutorial archives do **not** agree with each
other, so this is measured rather than assumed — see ``data.H_SCORE_ASCENDING``.

For every pair it labels each correspondence as a ground-truth inlier or not
(reprojection under H_gt for H, Sampson distance to the GT fundamental for F)
and computes the AUC of "score ranks inliers before outliers". AUC > 0.5 means
lower scores are better; < 0.5 means higher scores are better.

    python check_scores.py
"""
import numpy as np
import cv2

import data

TH_PX = 3.0


def _auc_lower_is_better(scores, is_inlier):
    """P(random inlier scores below a random outlier), by rank."""
    order = np.argsort(scores, kind="stable")
    rank = np.empty(len(scores))
    rank[order] = np.arange(len(scores))
    n1, n0 = is_inlier.sum(), (~is_inlier).sum()
    if n1 == 0 or n0 == 0:
        return None
    return 1 - (rank[is_inlier].sum() - n1 * (n1 - 1) / 2) / (n1 * n0)


def _h_inliers(pair, th=TH_PX):
    src = pair["pts1"].reshape(-1, 1, 2).astype(np.float32)
    proj = cv2.perspectiveTransform(
        src, pair["H_gt"].astype(np.float32)).reshape(-1, 2)
    return np.sqrt(((proj - pair["pts2"]) ** 2).sum(1)) <= th


def _gt_fundamental(pair):
    R = pair["R2"] @ pair["R1"].T
    t = (pair["T2"] - R @ pair["T1"]).ravel()
    tx = np.array([[0, -t[2], t[1]], [t[2], 0, -t[0]], [-t[1], t[0], 0]])
    return np.linalg.inv(pair["K2"]).T @ (tx @ R) @ np.linalg.inv(pair["K1"])


def _f_inliers(pair, th=TH_PX):
    F = _gt_fundamental(pair)
    x1 = np.hstack([pair["pts1"], np.ones((len(pair["pts1"]), 1))])
    x2 = np.hstack([pair["pts2"], np.ones((len(pair["pts2"]), 1))])
    Fx1, Ftx2 = (F @ x1.T).T, (F.T @ x2.T).T
    num = np.sum(x2 * Fx1, 1) ** 2
    den = Fx1[:, 0] ** 2 + Fx1[:, 1] ** 2 + Ftx2[:, 0] ** 2 + Ftx2[:, 1] ** 2
    return num / np.maximum(den, 1e-12) <= th ** 2


def summarise(name, pairs, inlier_fn, ascending, limit=0):
    """Report the orientation of the *raw* stored scores, and check that what
    the loaders hand out is normalised to lower-is-better."""
    raw, norm = [], []
    for i, pair in enumerate(pairs):
        if limit and i >= limit:
            break
        inl = inlier_fn(pair)
        auc = _auc_lower_is_better(pair["scores"], inl)
        if auc is None:
            continue
        norm.append(auc)
        # The loaders negate reversed datasets; negating again recovers the
        # raw stored score, since negation is its own inverse.
        raw.append(_auc_lower_is_better(
            data._oriented(pair["scores"], ascending), inl))
    r, n = np.array(raw), np.array(norm)
    verdict = "lower is better" if r.mean() > 0.5 else "HIGHER is better"
    print(f"{name:24s} pairs={len(r):4d}  raw AUC={r.mean():.3f}  "
          f"agreeing={np.mean(r > 0.5):.0%}  -> stored {verdict:16s}"
          f"| after loader: {n.mean():.3f} "
          f"{'ok' if n.mean() > 0.5 else 'STILL REVERSED'}")
    return r.mean() > 0.5


def main():
    print("AUC of 'low score ranks GT inliers first'. >0.5: lower is better.\n")
    results = {}
    for ds in data.H_DATASETS:
        results[ds] = summarise(f"H {ds}", data.iter_pairs_h(ds, "test"),
                                _h_inliers, data.H_SCORE_ASCENDING[ds])
    results["F"] = summarise("F " + data.F_SCENES[0], data.iter_pairs_f(),
                             _f_inliers, data.F_SCORE_ASCENDING, limit=200)

    print("\nconfigured in data.py:")
    for ds in data.H_DATASETS:
        ok = "ok" if data.H_SCORE_ASCENDING[ds] == results[ds] else "MISMATCH"
        print(f"  H_SCORE_ASCENDING[{ds!r}] = {data.H_SCORE_ASCENDING[ds]}"
              f"  ({ok})")
    print(f"  F scores are used ascending = {data.F_SCORE_ASCENDING} "
          f"({'ok' if data.F_SCORE_ASCENDING == results['F'] else 'MISMATCH'})")


if __name__ == "__main__":
    main()
