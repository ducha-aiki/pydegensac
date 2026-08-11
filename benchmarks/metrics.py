"""Accuracy metrics, ported from the reference evaluations.

F: the IMC pose metric — recover relative pose from the estimated fundamental
matrix plus ground-truth intrinsics, take max(rotation error, translation
direction error) in degrees, and average the fraction below each threshold
1..10 deg. Ported from imc2021-simple ``imc2021/metrics.py``
(``pose_error`` + ``maa_imc``), which matches ``qt_auc_10`` of the
image-matching-benchmark.

H: the CVPR-2020 RANSAC tutorial metric — reproject the image-1 pixel grid
under the ground-truth and the estimated homography, average the Euclidean
residual over the jointly visible area, and average the fraction below each of
10 log-spaced thresholds 1..20 px. Ported from the tutorial's ``metrics.py``
(``get_visible_part_mean_absolute_reprojection_error`` + ``calc_mAA``).
"""
import cv2
import numpy as np

#: IMC angular thresholds, degrees.
F_THRESHOLDS = list(range(1, 11))
#: Tutorial reprojection thresholds, pixels.
H_THRESHOLDS = np.logspace(np.log2(1.0), np.log2(20.0), 10, base=2.0)

#: Error assigned to a pair where the estimator returned no usable model.
#: Larger than any threshold, so mAA counts it as a miss at every level.
FAIL_ERR = np.inf


def pose_error(F, kp1, kp2, K1, K2, R1, R2, T1, T2):
    """max(R_err, t_err) in degrees for an estimated fundamental matrix.

    ``kp1``/``kp2`` are the inlier correspondences in pixel coordinates.
    Returns ``FAIL_ERR`` on degenerate input or a failed decomposition.
    """
    if len(kp1) < 5:
        return FAIL_ERR

    R_gt = R2 @ R1.T
    t_gt = (T2 - R_gt @ T1).flatten()
    if np.linalg.norm(t_gt) < 1e-9:
        return FAIL_ERR

    E = K2.T @ F @ K1
    norm = lambda kp, K: (kp - K[:2, 2]) / np.array([K[0, 0], K[1, 1]])  # noqa: E731
    try:
        _, R, t, _ = cv2.recoverPose(E, norm(kp1, K1), norm(kp2, K2))
    except cv2.error:
        return FAIL_ERR

    cos_R = np.clip((np.trace(R @ R_gt.T) - 1) / 2, -1, 1)
    R_err = np.degrees(np.arccos(cos_R))

    eps = 1e-15
    t = t.flatten()
    cos_t = np.abs(np.dot(t / (np.linalg.norm(t) + eps),
                          t_gt / (np.linalg.norm(t_gt) + eps)))
    t_err = np.degrees(np.arccos(np.clip(cos_t, 0, 1)))

    return float(max(R_err, t_err))


def reprojection_error(shape1, shape2, H_gt, H):
    """Mean reprojection error over the jointly visible part of image 1.

    ``shape1``/``shape2`` are (h, w). Returns ``FAIL_ERR`` when the estimate is
    unusable or the images share no visible area.
    """
    if H is None or not np.all(np.isfinite(H)):
        return FAIL_ERR
    h, w = shape1
    try:
        H_gt_inv = np.linalg.inv(H_gt)
    except np.linalg.LinAlgError:
        return FAIL_ERR
    mask1 = np.ones((h, w), np.float32)
    mask1in2 = cv2.warpPerspective(mask1, H_gt, (shape2[1], shape2[0]))
    visible = cv2.warpPerspective(mask1in2, H_gt_inv, (w, h)) > 0
    if not visible.any():
        return FAIL_ERR

    xg, yg = np.meshgrid(np.arange(w), np.arange(h))
    coords = np.stack([xg, yg], -1).reshape(-1, 1, 2).astype(np.float32)
    gt = cv2.perspectiveTransform(coords, H_gt.astype(np.float32)).squeeze(1)
    est = cv2.perspectiveTransform(coords, H.astype(np.float32)).squeeze(1)
    err = np.sqrt(((gt - est) ** 2).sum(1)).reshape(h, w) * visible
    if not np.isfinite(err).all():
        return FAIL_ERR
    return float(err.sum() / visible.sum())


def maa(errors, thresholds, strict=False):
    """Mean over thresholds of the fraction of pairs below each threshold.

    ``strict`` picks ``<`` over ``<=`` — the two references differ on this and
    each is followed exactly, though with continuous errors it never matters.
    """
    errors = np.asarray(errors, float)
    if errors.size == 0:
        return 0.0
    cmp = (lambda t: errors < t) if strict else (lambda t: errors <= t)
    return float(np.mean([cmp(t).mean() for t in thresholds]))


def maa_f(errors):
    """IMC mAA over 1..10 deg (imc2021-simple ``maa_imc``: strict)."""
    return maa(errors, F_THRESHOLDS, strict=True)


def maa_h(errors):
    """Tutorial mAA over 1..20 px (tutorial ``calc_mAA``: inclusive)."""
    return maa(errors, H_THRESHOLDS, strict=False)
