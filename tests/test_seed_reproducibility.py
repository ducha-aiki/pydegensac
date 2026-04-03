import numpy as np

import pydegensac


def _make_homography_data():
    rng = np.random.default_rng(0)
    pts1_in = rng.uniform(-50.0, 50.0, size=(48, 2))
    H_true = np.array(
        [
            [1.02, 0.04, 8.0],
            [-0.03, 0.97, -5.0],
            [0.0006, -0.0004, 1.0],
        ],
        dtype=np.float64,
    )

    homog = np.c_[pts1_in, np.ones(len(pts1_in))]
    proj = homog @ H_true.T
    pts2_in = proj[:, :2] / proj[:, 2:]

    pts1_out = rng.uniform(-50.0, 50.0, size=(24, 2))
    pts2_out = rng.uniform(-50.0, 50.0, size=(24, 2))

    pts1 = np.vstack([pts1_in, pts1_out])
    pts2 = np.vstack([pts2_in, pts2_out])
    return pts1, pts2


def _make_fundamental_data():
    rng = np.random.default_rng(1)
    pts3d = rng.uniform([-1.0, -1.0, 4.0], [1.0, 1.0, 8.0], size=(64, 3))

    p1 = np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
        ],
        dtype=np.float64,
    )
    p2 = np.array(
        [
            [0.999, 0.0, 0.045, 0.35],
            [0.002, 1.0, 0.0, -0.1],
            [-0.045, 0.0, 0.999, 0.2],
        ],
        dtype=np.float64,
    )

    homog = np.c_[pts3d, np.ones(len(pts3d))]
    x1 = homog @ p1.T
    x2 = homog @ p2.T
    pts1_in = x1[:, :2] / x1[:, 2:]
    pts2_in = x2[:, :2] / x2[:, 2:]

    pts1_out = rng.uniform(-0.5, 0.5, size=(24, 2))
    pts2_out = rng.uniform(-0.5, 0.5, size=(24, 2))

    pts1 = np.vstack([pts1_in, pts1_out])
    pts2 = np.vstack([pts2_in, pts2_out])
    return pts1, pts2


def test_seed_reproducibility_homography():
    pts1, pts2 = _make_homography_data()

    h1, mask1 = pydegensac.findHomography(
        pts1, pts2, px_th=0.25, conf=0.999, max_iters=5000, seed=1234
    )
    h2, mask2 = pydegensac.findHomography(
        pts1, pts2, px_th=0.25, conf=0.999, max_iters=5000, seed=1234
    )

    assert np.array_equal(np.asarray(mask1), np.asarray(mask2))
    assert np.allclose(h1, h2, rtol=0.0, atol=1e-12)


def test_seed_reproducibility_fundamental():
    pts1, pts2 = _make_fundamental_data()

    f1, mask1 = pydegensac.findFundamentalMatrix(
        pts1,
        pts2,
        px_th=1e-4,
        conf=0.999,
        max_iters=20000,
        enable_degeneracy_check=True,
        seed=5678,
    )
    f2, mask2 = pydegensac.findFundamentalMatrix(
        pts1,
        pts2,
        px_th=1e-4,
        conf=0.999,
        max_iters=20000,
        enable_degeneracy_check=True,
        seed=5678,
    )

    assert np.array_equal(np.asarray(mask1), np.asarray(mask2))
    assert np.allclose(f1, f2, rtol=0.0, atol=1e-9)


def test_seeded_homography_returns_valid_model():
    pts1, pts2 = _make_homography_data()

    h, mask = pydegensac.findHomography(
        pts1, pts2, px_th=0.25, conf=0.999, max_iters=50, seed=111
    )

    assert h.shape == (3, 3)
    assert len(mask) == len(pts1)
    assert np.isfinite(h).all()
    assert int(np.asarray(mask, dtype=bool).sum()) >= 48


def test_seeded_fundamental_returns_valid_model():
    pts1, pts2 = _make_fundamental_data()

    f, mask = pydegensac.findFundamentalMatrix(
        pts1,
        pts2,
        px_th=1e-4,
        conf=0.999,
        max_iters=200,
        enable_degeneracy_check=True,
        seed=111,
    )

    assert f.shape == (3, 3)
    assert len(mask) == len(pts1)
    assert np.isfinite(f).all()
    assert int(np.asarray(mask, dtype=bool).sum()) >= 64
