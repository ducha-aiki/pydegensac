import numpy as np
import pytest

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


def test_public_api_seed_makes_results_deterministic():
    pts1_h, pts2_h, _ = _make_homography_points()
    h1, mask1 = pydegensac.findHomography(pts1_h, pts2_h, max_iters=100, seed=123)
    h2, mask2 = pydegensac.findHomography(pts1_h, pts2_h, max_iters=100, seed=123)
    np.testing.assert_array_equal(mask1, mask2)
    np.testing.assert_allclose(h1, h2)

    pts1_f, pts2_f = _make_fundamental_points()
    f1, inliers1 = pydegensac.findFundamentalMatrix(
        pts1_f,
        pts2_f,
        max_iters=1000,
        enable_degeneracy_check=False,
        seed=123,
    )
    f2, inliers2 = pydegensac.findFundamentalMatrix(
        pts1_f,
        pts2_f,
        max_iters=1000,
        enable_degeneracy_check=False,
        seed=123,
    )
    np.testing.assert_array_equal(inliers1, inliers2)
    np.testing.assert_allclose(f1, f2)


def test_concurrent_threads_deterministic_with_degeneracy_check():
    from concurrent.futures import ThreadPoolExecutor

    pts1_f, pts2_f = _make_fundamental_points()

    # Fundamental matrix with degeneracy check enabled
    f_serial, m_serial = pydegensac.findFundamentalMatrix(
        pts1_f, pts2_f, px_th=1e-4, enable_degeneracy_check=True, max_iters=1000, seed=42
    )

    with ThreadPoolExecutor(max_workers=8) as ex:
        threaded_f = list(
            ex.map(
                lambda _: pydegensac.findFundamentalMatrix(
                    pts1_f, pts2_f, px_th=1e-4, enable_degeneracy_check=True, max_iters=1000, seed=42
                ),
                range(20),
            )
        )

    for f_par, m_par in threaded_f:
        np.testing.assert_array_equal(m_serial, m_par)
        np.testing.assert_allclose(f_serial, f_par)

    # Homography matrix under concurrent threads
    pts1_h, pts2_h, _ = _make_homography_points()
    h_serial, m_h_serial = pydegensac.findHomography(
        pts1_h, pts2_h, px_th=1e-5, max_iters=1000, seed=42
    )

    with ThreadPoolExecutor(max_workers=8) as ex:
        threaded_h = list(
            ex.map(
                lambda _: pydegensac.findHomography(
                    pts1_h, pts2_h, px_th=1e-5, max_iters=1000, seed=42
                ),
                range(20),
            )
        )

    for h_par, m_par in threaded_h:
        np.testing.assert_array_equal(m_h_serial, m_par)
        np.testing.assert_allclose(h_serial, h_par)


def test_homography_laf_consistency():
    rng = np.random.RandomState(42)
    pts1 = rng.uniform(50, 500, size=(40, 2))
    theta = 0.05
    r_mat = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    t_vec = np.array([20.0, 30.0])

    pts2 = (pts1 @ r_mat.T) + t_vec

    laf1 = np.zeros((len(pts1), 6), dtype=np.float64)
    laf1[:, :2] = pts1
    laf1[:, 2] = 1.0
    laf1[:, 5] = 1.0

    laf2 = np.zeros((len(pts2), 6), dtype=np.float64)
    laf2[:, :2] = pts2
    laf2[:, 2] = r_mat[0, 0]
    laf2[:, 3] = r_mat[0, 1]
    laf2[:, 4] = r_mat[1, 0]
    laf2[:, 5] = r_mat[1, 1]

    # Invert affine orientation for the last 10 points
    laf2[30:, 2:6] = -laf2[30:, 2:6]

    # Without LAF check, all 40 points have consistent positions
    _, mask_no_laf = pydegensac.findHomography(
        laf1, laf2, px_th=1.0, laf_consistensy_coef=-1.0, seed=42
    )
    assert np.count_nonzero(mask_no_laf) == 40

    # With LAF consistency check, the 10 affine-inconsistent points are rejected
    h_laf, mask_laf = pydegensac.findHomography(
        laf1, laf2, px_th=1.0, laf_consistensy_coef=1.0, seed=42
    )
    assert np.count_nonzero(mask_laf[:30]) == 30
    assert np.count_nonzero(mask_laf[30:]) == 0
    assert h_laf.shape == (3, 3)
    assert np.isfinite(h_laf).all()


@pytest.mark.parametrize("estimator", ["homography", "fundamental"])
@pytest.mark.parametrize("sym_check,laf_coef", [(False, 0.0), (True, 0.0),
                                             (False, 1.0), (True, 1.0)])
def test_native_mask_owns_storage_across_calls(estimator, sym_check, laf_coef):
    # Odd-sized inputs also exercise alignment between double and int scratch.
    rng = np.random.RandomState(19)
    pts1 = rng.uniform(0, 100, (65, 2))
    pts2 = pts1 + [3.0, -2.0]
    if estimator == "fundamental":
        xyz = rng.uniform([-2, -2, 4], [2, 2, 10], (65, 3))
        pts1 = 800 * xyz[:, :2] / xyz[:, 2:]
        xyz2 = xyz + [0.5, 0.1, 0.05]
        pts2 = 800 * xyz2[:, :2] / xyz2[:, 2:]

    if laf_coef:
        frame = np.tile([1.0, 0.0, 0.0, 1.0], (len(pts1), 1))
        pts1 = np.column_stack((pts1, frame))
        pts2 = np.column_stack((pts2, frame))
    if estimator == "homography":
        call = lambda: pydegensac.findHomography_(
            pts1, pts2, 1.0, 0.999, 10, 0, sym_check, laf_coef, 42)
    else:
        call = lambda: pydegensac.findFundamentalMatrix_(
            pts1, pts2, 1.0, 0.999, 10, 0, sym_check, laf_coef, True, 42)
    model, mask = call()
    assert mask.any()
    saved_model, saved_mask = model.copy(), mask.copy()
    for _ in range(3):
        next_model, next_mask = call()
        assert next_mask.dtype == np.bool_
        assert next_mask.shape == (len(pts1),)
        np.testing.assert_array_equal(next_model, saved_model)
        np.testing.assert_array_equal(next_mask, saved_mask)
        next_mask[:] = False
        next_model[:] = 0
        np.testing.assert_array_equal(model, saved_model)
        np.testing.assert_array_equal(mask, saved_mask)


def test_degenerate_native_calls_return_empty_masks():
    points = np.zeros((8, 2), dtype=np.float64)

    h, h_mask = pydegensac.findHomography_(
        points[:4], points[:4], 1.0, 0.999, 1, 0, False, 0.0, 42
    )
    f, f_mask = pydegensac.findFundamentalMatrix_(
        points, points, 1.0, 0.999, 1, 0, False, 0.0, False, 42
    )

    np.testing.assert_array_equal(h, np.zeros((3, 3)))
    np.testing.assert_array_equal(f, np.zeros((3, 3)))
    assert not h_mask.any()
    assert not f_mask.any()
