import numpy as np

import pydegensac


def _make_homography_points():
    pts1 = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [0.2, 0.7],
            [0.8, 0.3],
            [0.4, 0.6],
            [0.6, 0.2],
        ],
        dtype=np.float64,
    )
    h_true = np.array(
        [[1.2, 0.1, 0.5], [0.05, 0.9, -0.2], [0.001, 0.002, 1.0]],
        dtype=np.float64,
    )
    pts1_h = np.concatenate([pts1, np.ones((len(pts1), 1), dtype=np.float64)], axis=1)
    pts2_h = (h_true @ pts1_h.T).T
    pts2 = pts2_h[:, :2] / pts2_h[:, 2:3]
    return pts1, pts2, h_true


def _make_fundamental_points():
    pts3d = np.array(
        [
            [0.1, 0.2, 2.0],
            [0.5, -0.1, 3.0],
            [-0.3, 0.4, 4.0],
            [0.2, 0.3, 2.5],
            [0.7, 0.5, 5.0],
            [-0.4, -0.2, 3.5],
            [0.0, 0.1, 2.2],
            [0.3, -0.4, 4.5],
            [0.6, 0.2, 3.2],
            [-0.2, 0.6, 5.5],
        ],
        dtype=np.float64,
    )
    pts1 = pts3d[:, :2] / pts3d[:, 2:3]
    pts2 = (pts3d[:, :2] + np.array([1.0, 0.0])) / pts3d[:, 2:3]
    return pts1, pts2


def _project_homography(h, pts):
    pts_h = np.concatenate([pts, np.ones((len(pts), 1), dtype=np.float64)], axis=1)
    out_h = (h @ pts_h.T).T
    return out_h[:, :2] / out_h[:, 2:3]


def test_repeated_homography_calls_return_valid_models():
    pts1, pts2, _ = _make_homography_points()

    for _ in range(5):
        h, mask = pydegensac.findHomography(pts1, pts2, px_th=1e-6, conf=0.999, max_iters=100)
        assert h.shape == (3, 3)
        assert np.isfinite(h).all()
        assert mask.dtype == np.bool_
        assert mask.shape == (len(pts1),)
        assert np.count_nonzero(mask) >= 4

        h_low, mask_low = pydegensac.findHomography_(
            pts1,
            pts2,
            1e-6,
            0.999,
            100,
            0,
            True,
            0.0,
            123,
        )
        assert h_low.shape == (3, 3)
        assert np.isfinite(h_low).all()
        assert mask_low.shape == (len(pts1),)
        assert np.count_nonzero(mask_low) >= 4


def test_repeated_fundamental_calls_return_epipolar_consistency():
    pts1, pts2 = _make_fundamental_points()
    pts1_h = np.concatenate([pts1, np.ones((len(pts1), 1), dtype=np.float64)], axis=1)
    pts2_h = np.concatenate([pts2, np.ones((len(pts2), 1), dtype=np.float64)], axis=1)

    for _ in range(5):
        f, mask = pydegensac.findFundamentalMatrix(
            pts1,
            pts2,
            px_th=1e-6,
            conf=0.999,
            max_iters=1000,
            enable_degeneracy_check=False,
        )
        residuals = np.abs(np.sum((pts2_h @ f) * pts1_h, axis=1))
        assert mask.dtype == np.bool_
        assert mask.shape == (len(pts1),)
        assert mask.all()
        assert np.max(residuals) < 1e-6

        f_low, mask_low = pydegensac.findFundamentalMatrix_(
            pts1,
            pts2,
            1e-6,
            0.999,
            1000,
            0,
            True,
            0.0,
            False,
            123,
        )
        residuals_low = np.abs(np.sum((pts2_h @ f_low) * pts1_h, axis=1))
        assert f_low.shape == (3, 3)
        assert np.isfinite(f_low).all()
        assert mask_low.shape == (len(pts1),)
        assert mask_low.all()
        assert np.max(residuals_low) < 1e-6


def test_public_api_forwards_seed_to_low_level_bindings(monkeypatch):
    seen = {}

    def fake_find_h(pts1, pts2, px_th, conf, max_iters, error_type, symmetric_error_check, laf_coef, seed):
        seen["h_seed"] = seed
        return np.eye(3, dtype=np.float64), np.ones(len(pts1), dtype=bool)

    def fake_find_f(
        pts1,
        pts2,
        px_th,
        conf,
        max_iters,
        error_type,
        symmetric_error_check,
        laf_coef,
        enable_degeneracy_check,
        seed,
    ):
        seen["f_seed"] = seed
        return np.eye(3, dtype=np.float64), np.ones(len(pts1), dtype=bool)

    monkeypatch.setattr(pydegensac, "findHomography_", fake_find_h)
    monkeypatch.setattr(pydegensac, "findFundamentalMatrix_", fake_find_f)

    pts1_h, pts2_h, _ = _make_homography_points()
    pydegensac.findHomography(pts1_h, pts2_h, max_iters=100, seed=123)
    assert seen["h_seed"] == 123

    pts1_f, pts2_f = _make_fundamental_points()
    pydegensac.findFundamentalMatrix(
        pts1_f,
        pts2_f,
        max_iters=1000,
        enable_degeneracy_check=False,
        seed=456,
    )
    assert seen["f_seed"] == 456
