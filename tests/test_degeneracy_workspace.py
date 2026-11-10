"""Repeated planar-scene recovery must not retain another call's scratch."""
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest

import pydegensac


def _scene(n, plane_fraction):
    rng = np.random.RandomState(n)
    xyz = rng.uniform([-2, -2, 2.5], [2, 2, 4], (n, 3))
    xyz[:int(n * plane_fraction), 2] = 6
    angle = 0.08
    rotation = np.array([[np.cos(angle), 0, np.sin(angle)],
                         [0, 1, 0], [-np.sin(angle), 0, np.cos(angle)]])
    xyz2 = xyz @ rotation.T + [0.8, 0.1, 0.05]
    pts1 = 800 * xyz[:, :2] / xyz[:, 2:]
    pts2 = 800 * xyz2[:, :2] / xyz2[:, 2:]
    pts1 += rng.normal(0, 0.1, pts1.shape)
    pts2 += rng.normal(0, 0.1, pts2.shape)
    return pts1, pts2


@pytest.mark.parametrize("error_type", ["sampson", "symm_epipolar"])
@pytest.mark.parametrize("plane_fraction", [0.6, 0.9, 1.0])
def test_planar_recovery_reuses_scratch_safely(plane_fraction, error_type):
    """Recovery remains deterministic with and without off-plane support."""
    cases = [(*_scene(n, plane_fraction), seed)
             for n in (65, 257) for seed in (42, 2024)]

    def estimate(case):
        pts1, pts2, seed = case
        return pydegensac.findFundamentalMatrix(
            pts1, pts2, px_th=0.75, max_iters=100, seed=seed,
            error_type=error_type, enable_degeneracy_check=True)

    expected = [estimate(case) for case in cases]
    for (model, mask), (pts1, _, _) in zip(expected, cases):
        assert np.isfinite(model).all()
        assert mask.sum() >= 0.8 * len(pts1)
    with ThreadPoolExecutor(max_workers=4) as executor:
        actual = list(executor.map(estimate, cases * 3))
    for (model, mask), (ref_model, ref_mask) in zip(actual, expected * 3):
        np.testing.assert_array_equal(model, ref_model)
        np.testing.assert_array_equal(mask, ref_mask)
