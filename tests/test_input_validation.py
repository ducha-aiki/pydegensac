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
