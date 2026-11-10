"""Exact regression tests for the deferred Sampson LO weights."""
import ctypes
import os

import numpy as np
import pytest

import pydegensac


_DOUBLE_PTR = ctypes.POINTER(ctypes.c_double)
_INT_PTR = ctypes.POINTER(ctypes.c_int)
_BYTE_PTR = ctypes.POINTER(ctypes.c_ubyte)
_ERROR_CALLBACK = ctypes.CFUNCTYPE(None, _DOUBLE_PTR, _DOUBLE_PTR, _DOUBLE_PTR,
                                   _DOUBLE_PTR, ctypes.c_int)


@pytest.fixture(scope="module")
def native_helpers():
    """Load C helpers from the extension when its symbols are exported."""
    try:
        library = ctypes.CDLL(pydegensac.pydegensac.__file__)
        residuals = library.exFDs
        denominators = library.exFDs_den
        selected_weights = library.exFDs_weights
    except (AttributeError, OSError):
        if os.name == "nt":
            pytest.skip("Windows extension does not export LO-weight helpers")
        raise

    for function in (residuals, denominators):
        function.argtypes = [_DOUBLE_PTR, _DOUBLE_PTR, _DOUBLE_PTR, _DOUBLE_PTR,
                             ctypes.c_int]
        function.restype = None
    selected_weights.argtypes = [_DOUBLE_PTR, _INT_PTR, ctypes.c_int]
    selected_weights.restype = None
    return residuals, denominators, selected_weights


def _ptr(array, pointer_type):
    return array.ctypes.data_as(pointer_type)


def _call_error(function, points, matrix, residuals, weights):
    function(_ptr(points, _DOUBLE_PTR), _ptr(matrix, _DOUBLE_PTR),
             _ptr(residuals, _DOUBLE_PTR), _ptr(weights, _DOUBLE_PTR), len(points))


def _bytes(array):
    return np.ascontiguousarray(array).view(np.uint8).tobytes()


def _points(length, kind):
    if kind == "zero":
        return np.zeros((length, 6), dtype=np.float64)
    if kind == "signedzero":
        points = np.zeros((length, 6), dtype=np.float64)
        points.flat[1::2] = -0.0
        return points
    if kind == "nonfinite":
        points = np.zeros((length, 6), dtype=np.float64)
        if length:
            points[0, 0] = np.nan
        if length > 1:
            points[1, 3] = np.inf
        return points
    rng = np.random.RandomState(91 + length)
    points = rng.uniform(-1000, 1000, (length, 6)).astype(np.float64)
    # The two homogeneous coordinates are not used by Sampson residuals, but
    # keep fixtures representative of converted (x, y, 1, x', y', 1) data.
    points[:, 2] = 1.0
    points[:, 5] = 1.0
    if kind == "extreme":
        points *= 1e150
    return points


def _matrix(kind):
    if kind == "zero":
        return np.zeros(9, dtype=np.float64)
    if kind == "signedzero":
        matrix = np.zeros(9, dtype=np.float64)
        matrix[1::2] = -0.0
        return matrix
    if kind == "nonfinite":
        matrix = np.arange(1, 10, dtype=np.float64)
        matrix[3] = np.nan
        return matrix
    matrix = np.array([0.001, -0.002, 0.3, 0.004, 0.005, -0.2,
                       0.8, -0.6, 1.0], dtype=np.float64)
    if kind == "extreme":
        matrix *= 1e150
    return matrix


@pytest.mark.parametrize("length", [0, 1, 7, 8, 9, 16, 17, 65, 257, 1024])
@pytest.mark.parametrize("kind", ["ordinary", "zero", "signedzero", "extreme", "nonfinite"])
def test_deferred_sampson_weights_match_full_error(native_helpers, length, kind):
    full_error, denominator_error, convert_weights = native_helpers
    points = _points(length, kind)
    matrix = _matrix(kind)
    full_residuals = np.full(length, -1.0, dtype=np.float64)
    full_weights = np.full(length, -2.0, dtype=np.float64)
    deferred_residuals = np.full(length, -3.0, dtype=np.float64)
    deferred_weights = np.full(length, -4.0, dtype=np.float64)

    _call_error(full_error, points, matrix, full_residuals, full_weights)
    _call_error(denominator_error, points, matrix, deferred_residuals, deferred_weights)
    assert _bytes(deferred_residuals) == _bytes(full_residuals)

    choices = [np.empty(0, dtype=np.int32)]
    if length >= 8:
        choices.append(np.arange(0, 8, dtype=np.int32))
        # randsubset returns a unique tail whose order is not necessarily
        # ascending; weights must be independent of that order.
        choices.append(np.arange(length - 1, length - 9, -1, dtype=np.int32))
    choices.append(np.arange(length, dtype=np.int32))
    for selected in choices:
        candidate = deferred_weights.copy()
        before = candidate.copy()
        convert_weights(_ptr(candidate, _DOUBLE_PTR), _ptr(selected, _INT_PTR), len(selected))
        expected = full_weights.copy()
        unselected = np.ones(length, dtype=bool)
        unselected[selected] = False
        assert _bytes(candidate[selected]) == _bytes(expected[selected])
        assert _bytes(candidate[unselected]) == _bytes(before[unselected])

    # Reusing the scratch after a partially converted fit must restore all
    # denominators exactly for the next LO refinement.
    if length:
        selected = np.arange(min(8, length), dtype=np.int32)
        convert_weights(_ptr(deferred_weights, _DOUBLE_PTR), _ptr(selected, _INT_PTR), len(selected))
        _call_error(denominator_error, points, matrix, deferred_residuals, deferred_weights)
        fresh_denominators = np.empty(length, dtype=np.float64)
        fresh_residuals = np.empty(length, dtype=np.float64)
        _call_error(denominator_error, points, matrix, fresh_residuals, fresh_denominators)
        assert _bytes(deferred_residuals) == _bytes(fresh_residuals)
        assert _bytes(deferred_weights) == _bytes(fresh_denominators)


def _as_void_pointer(value):
    if isinstance(value, np.ndarray):
        return ctypes.c_void_p(value.ctypes.data)
    return ctypes.cast(value, ctypes.c_void_p)


@pytest.mark.parametrize("inlier_limit", [0, 8, 49, 1000000])
def test_selected_weight_path_matches_custom_error_callback(native_helpers, inlier_limit):
    """The optimized built-in callback must match the legacy callback path."""
    full_error, _, _ = native_helpers
    library = ctypes.CDLL(pydegensac.pydegensac.__file__)
    try:
        ransac = library.exp_ransacFcustomLAF
        fds = library.FDs
        fdsidx = library.FDsidx
        fds_soa = library.FDs_soa
        library.degensac_rand.restype = ctypes.c_int
        library.degensac_random.restype = ctypes.c_long
    except AttributeError:
        pytest.skip("the native extension does not export the F RANSAC callbacks")

    # Keeping pointers opaque avoids platform-specific ctypes function-pointer
    # declarations while still exercising C's callback dispatch exactly.
    ransac.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
                       ctypes.c_int, ctypes.c_double, ctypes.c_double,
                       ctypes.c_double, ctypes.c_int, ctypes.c_void_p,
                       ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int,
                       ctypes.c_uint, ctypes.c_void_p, ctypes.c_void_p,
                       ctypes.c_void_p, ctypes.c_void_p,
                       ctypes.c_double, ctypes.c_int, ctypes.c_int,
                       ctypes.c_void_p]
    ransac.restype = ctypes.c_int

    rng = np.random.RandomState(13)
    xyz = rng.uniform([-2, -2, 4], [2, 2, 10], (65, 3))
    points1 = 800 * xyz[:, :2] / xyz[:, 2:]
    points2 = 800 * (xyz[:, :2] + [0.5, 0.1]) / (xyz[:, 2:] + 0.05)
    points1 += rng.normal(0, 0.05, points1.shape)
    points2 += rng.normal(0, 0.05, points2.shape)
    points2[:12] = rng.uniform(-400, 400, (12, 2))
    u = np.ascontiguousarray(np.column_stack((points1, np.ones(65),
                                               points2, np.ones(65))), dtype=np.float64)

    def run(error_callback):
        matrix = np.zeros(9, dtype=np.float64)
        mask = np.zeros(len(u), dtype=np.uint8)
        result = ransac(
            _as_void_pointer(u), None, None, len(u), 0.75 ** 2, 0.0, 0.999,
            100, _as_void_pointer(matrix), _as_void_pointer(mask), None,
            1, inlier_limit, None, _as_void_pointer(error_callback),
            _as_void_pointer(fds), _as_void_pointer(fdsidx), 0.0, 0, 42,
            _as_void_pointer(fds_soa),
        )
        return result, matrix, mask, library.degensac_rand(), library.degensac_random()

    calls = 0

    @_ERROR_CALLBACK
    def legacy_callback(u_ptr, f_ptr, residual_ptr, weight_ptr, length):
        nonlocal calls
        calls += 1
        full_error(u_ptr, f_ptr, residual_ptr, weight_ptr, length)

    built_in = run(full_error)
    legacy = run(legacy_callback)
    assert built_in[0] == legacy[0]
    assert _bytes(built_in[1]) == _bytes(legacy[1])
    assert _bytes(built_in[2]) == _bytes(legacy[2])
    assert built_in[3:] == legacy[3:]
    assert calls > 0
    assert np.any(built_in[1])
    assert np.count_nonzero(built_in[2]) >= 8
