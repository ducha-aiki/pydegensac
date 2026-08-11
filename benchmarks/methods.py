"""Uniform estimator wrappers.

Every wrapper has the signature ``fn(pts1, pts2, th, scores) -> (M, mask)``,
where ``M`` is None when the estimator produced nothing usable and ``mask`` is
a boolean inlier mask over the *input* order. ``scores`` are per-match
qualities, lower = better (SNN ratio); only the PROSAC variants read them.

The poselib PROSAC wrapper is ported from imc2021-simple
``imc2021/ransac.py``: poselib does not sort internally, so points go in
sorted best-first and the returned mask is un-permuted afterwards.

Confidence follows each reference evaluation: 0.9999 for F (imc2021-simple),
0.999 for H (ds-sac / the tutorial's own scripts).
"""
import cv2
import numpy as np

F_CONF = 0.9999
H_CONF = 0.999

#: Iteration budgets for the time-mAA sweeps.
F_BUDGETS = (125, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000)
H_BUDGETS = (10, 25, 100, 400, 1600, 6400, 25000)

METHOD_NAMES = ("cv2-ransac", "cv2-magsac", "poselib", "poselib-prosac",
                "pydegensac")


def _fail(n):
    return None, np.zeros(n, bool)


def _prosac_order(scores, n, use_prosac):
    """Permutation putting the best matches first, or None for uniform."""
    if not use_prosac or scores is None:
        return None
    return np.argsort(np.asarray(scores).ravel(), kind="stable")


def _unpermute(mask_sorted, order, n):
    mask = np.zeros(n, bool)
    mask[order] = mask_sorted
    return mask


# --------------------------------------------------------------------------
# Fundamental matrix
# --------------------------------------------------------------------------

def _cv2_f(flag):
    def run(pts1, pts2, th, scores=None, max_iters=F_BUDGETS[-1]):
        if len(pts1) < 8:
            return _fail(len(pts1))
        try:
            F, mask = cv2.findFundamentalMat(
                pts1, pts2, flag, th, confidence=F_CONF, maxIters=max_iters)
        except cv2.error:
            return _fail(len(pts1))
        if F is None or mask is None or F.shape[0] < 3:
            return _fail(len(pts1))
        # cv2 can return several stacked candidates; the first is the answer.
        return F[:3].astype(np.float64), mask.ravel().astype(bool)
    return run


def _poselib_f(use_prosac):
    def run(pts1, pts2, th, scores=None, max_iters=F_BUDGETS[-1]):
        import poselib
        n = len(pts1)
        if n < 8:
            return _fail(n)
        order = _prosac_order(scores, n, use_prosac)
        src = np.ascontiguousarray(pts1 if order is None else pts1[order])
        dst = np.ascontiguousarray(pts2 if order is None else pts2[order])
        try:
            F, info = poselib.estimate_fundamental(
                src, dst,
                {"max_epipolar_error": th,
                 "progressive_sampling": order is not None,
                 "max_iterations": max_iters,
                 "success_prob": F_CONF},
                {})
        except Exception:
            return _fail(n)
        if F is None:
            return _fail(n)
        mask = np.asarray(info["inliers"], bool)
        if order is not None:
            mask = _unpermute(mask, order, n)
        return np.asarray(F, np.float64), mask
    return run


def _pydegensac_f(pts1, pts2, th, scores=None, max_iters=F_BUDGETS[-1]):
    import pydegensac
    n = len(pts1)
    if n < 8:
        return _fail(n)
    try:
        F, mask = pydegensac.findFundamentalMatrix(
            pts1, pts2, th, conf=F_CONF, max_iters=max_iters)
    except Exception:
        return _fail(n)
    if F is None or not np.any(F):
        return _fail(n)
    return np.asarray(F, np.float64), np.asarray(mask, bool).ravel()


F_METHODS = {
    "cv2-ransac": _cv2_f(cv2.FM_RANSAC),
    "cv2-magsac": _cv2_f(cv2.USAC_MAGSAC),
    "poselib": _poselib_f(use_prosac=False),
    "poselib-prosac": _poselib_f(use_prosac=True),
    "pydegensac": _pydegensac_f,
}


# --------------------------------------------------------------------------
# Homography
# --------------------------------------------------------------------------

def _cv2_h(flag):
    def run(pts1, pts2, th, scores=None, max_iters=H_BUDGETS[-1]):
        if len(pts1) < 4:
            return _fail(len(pts1))
        try:
            H, mask = cv2.findHomography(
                pts1, pts2, flag, th, maxIters=max_iters, confidence=H_CONF)
        except cv2.error:
            return _fail(len(pts1))
        if H is None or mask is None:
            return _fail(len(pts1))
        return H.astype(np.float64), mask.ravel().astype(bool)
    return run


def _poselib_h(use_prosac):
    def run(pts1, pts2, th, scores=None, max_iters=H_BUDGETS[-1]):
        import poselib
        n = len(pts1)
        if n < 4:
            return _fail(n)
        order = _prosac_order(scores, n, use_prosac)
        src = np.ascontiguousarray(pts1 if order is None else pts1[order])
        dst = np.ascontiguousarray(pts2 if order is None else pts2[order])
        try:
            H, info = poselib.estimate_homography(
                src, dst,
                {"max_reproj_error": th,
                 "progressive_sampling": order is not None,
                 "max_iterations": max_iters,
                 "success_prob": H_CONF},
                {})
        except Exception:
            return _fail(n)
        if H is None:
            return _fail(n)
        mask = np.asarray(info["inliers"], bool)
        if order is not None:
            mask = _unpermute(mask, order, n)
        return np.asarray(H, np.float64), mask
    return run


def _pydegensac_h(pts1, pts2, th, scores=None, max_iters=H_BUDGETS[-1]):
    import pydegensac
    n = len(pts1)
    if n < 4:
        return _fail(n)
    try:
        H, mask = pydegensac.findHomography(
            pts1, pts2, th, conf=H_CONF, max_iters=max_iters)
    except Exception:
        return _fail(n)
    if H is None or not np.any(H):
        return _fail(n)
    return np.asarray(H, np.float64), np.asarray(mask, bool).ravel()


H_METHODS = {
    "cv2-ransac": _cv2_h(cv2.RANSAC),
    "cv2-magsac": _cv2_h(cv2.USAC_MAGSAC),
    "poselib": _poselib_h(use_prosac=False),
    "poselib-prosac": _poselib_h(use_prosac=True),
    "pydegensac": _pydegensac_h,
}


def registry(problem):
    return F_METHODS if problem == "f" else H_METHODS


def budgets(problem):
    return F_BUDGETS if problem == "f" else H_BUDGETS


def versions():
    """Backend versions, recorded with every result file.

    Both pydegensac builds report the same version string, so the install
    path is recorded too — that is what distinguishes the A/B arms.
    """
    import importlib
    import importlib.metadata as md
    out = {}
    for name in ("cv2", "poselib", "pydegensac", "numpy"):
        try:
            mod = importlib.import_module(name)
        except ImportError:
            out[name] = "missing"
            continue
        version = getattr(mod, "__version__", None)
        if version is None:
            try:
                version = md.version(name)
            except md.PackageNotFoundError:
                version = "unknown"
        out[name] = str(version)
    try:
        out["pydegensac_path"] = importlib.import_module("pydegensac").__file__
    except ImportError:
        pass
    return out
