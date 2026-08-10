import numpy as np
import pytest

import pydegensac
from pydegensac import utils


def test_failure_path_mask_is_bool_ndarray(monkeypatch):
    n = 8
    pts = np.random.RandomState(0).rand(n, 2)
    # force the no-model branch by stubbing the raw binding
    monkeypatch.setattr(utils.pydegensac, "findHomography_",
                        lambda *a, **k: (np.zeros((3, 3)), np.zeros(n)))
    H, mask = pydegensac.findHomography(pts, pts + 1)
    assert isinstance(mask, np.ndarray) and mask.dtype == np.bool_
    assert mask.shape == (n,) and not mask.any()
    monkeypatch.setattr(utils.pydegensac, "findFundamentalMatrix_",
                        lambda *a, **k: (np.zeros((3, 3)), np.zeros(n)))
    F, maskf = pydegensac.findFundamentalMatrix(pts, pts + 1)
    assert isinstance(maskf, np.ndarray) and maskf.dtype == np.bool_


def test_bad_error_type_raises_valueerror():
    pts = np.random.RandomState(0).rand(8, 2)
    with pytest.raises(ValueError, match="Error type"):
        pydegensac.findHomography(pts, pts, error_type="nope")
    with pytest.raises(ValueError, match="Error type"):
        pydegensac.findFundamentalMatrix(pts, pts, error_type="nope")


def test_fortran_order_gives_same_result_as_c_order():
    rs = np.random.RandomState(1)
    pts1 = rs.rand(60, 2) * 100
    H_true = np.array([[1.1, 0.05, 3.0], [0.02, 0.95, -2.0], [1e-4, 2e-4, 1.0]])
    p1h = np.concatenate([pts1, np.ones((60, 1))], axis=1)
    p2h = (H_true @ p1h.T).T
    pts2 = p2h[:, :2] / p2h[:, 2:3]
    Hc, mc = pydegensac.findHomography(pts1, pts2, 1.0, seed=7)
    Hf, mf = pydegensac.findHomography(np.asfortranarray(pts1),
                                       np.asfortranarray(pts2), 1.0, seed=7)
    np.testing.assert_array_equal(Hc, Hf)
    np.testing.assert_array_equal(np.asarray(mc), np.asarray(mf))


def test_1d_input_raises():
    with pytest.raises((ValueError, RuntimeError)):
        pydegensac.findHomography_(np.zeros(8), np.zeros(8),
                                   1.0, 0.999, 100, 0, True, 0.0, -1)


def test_mismatched_width_raises():
    with pytest.raises((ValueError, RuntimeError)):
        pydegensac.findHomography_(np.zeros((10, 6)), np.zeros((10, 2)),
                                   1.0, 0.999, 100, 0, True, 0.0, -1)


def test_laf_coef_with_2d_input_raises():
    pts = np.random.RandomState(2).rand(10, 2)
    with pytest.raises((ValueError, RuntimeError)):
        pydegensac.findHomography_(pts, pts, 1.0, 0.999, 100, 0, True, 0.5, -1)
    with pytest.raises((ValueError, RuntimeError)):
        pydegensac.findFundamentalMatrix_(pts, pts, 1.0, 0.999, 100, 0, True,
                                          0.5, True, -1)
