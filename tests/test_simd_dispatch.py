"""The AVX and scalar kernels must give bit-identical estimator output."""
import os
import platform
import subprocess
import sys
from pathlib import Path

import pytest

import pydegensac

DATA_DIR = Path(__file__).parent / "data"

# Runs in a fresh interpreter, because the backend is chosen once per process.
DIGEST = r"""
import hashlib, sys
from pathlib import Path
import numpy as np
import pydegensac

h = hashlib.sha256()
rng = np.random.default_rng(0)
for path in sorted(Path(sys.argv[1]).glob("golden_*.npz")):
    d = np.load(path)
    for seed in (1, 42, 2024):
        if path.name.startswith("golden_h"):
            for err in ("sampson", "symm_sq_max", "symm_max", "symm_sq_sum", "symm_sum"):
                M, mask = pydegensac.findHomography(d["pts1"], d["pts2"], 3.0, 0.99, 2000,
                                                    error_type=err, seed=seed)
                h.update(M.tobytes()); h.update(np.asarray(mask, bool).tobytes())
        else:
            for err in ("sampson", "symm_epipolar"):
                M, mask = pydegensac.findFundamentalMatrix(d["pts1"], d["pts2"], 0.5, 0.999, 5000,
                                                           error_type=err, seed=seed)
                h.update(M.tobytes()); h.update(np.asarray(mask, bool).tobytes())
# Lengths that are not multiples of four exercise the scalar tails.
for n in (8, 9, 10, 11, 37, 101):
    pts1 = rng.uniform(0, 640, (n, 2))
    pts2 = pts1 + rng.normal(0, 1, (n, 2))
    M, mask = pydegensac.findFundamentalMatrix(pts1, pts2, 1.0, 0.999, 500, seed=7)
    h.update(M.tobytes()); h.update(np.asarray(mask, bool).tobytes())
print(pydegensac.simd_backend_(), h.hexdigest())
"""


def _run(simd=None):
    env = dict(os.environ)
    env.pop("PYDEGENSAC_SIMD", None)
    if simd is not None:
        env["PYDEGENSAC_SIMD"] = simd
    out = subprocess.run([sys.executable, "-c", DIGEST, str(DATA_DIR)], env=env,
                         check=True, capture_output=True, text=True).stdout.split()
    return out[0], out[1]


def test_backend_name():
    assert pydegensac.simd_backend_() in ("avx", "scalar")


def test_env_forces_scalar():
    backend, _ = _run("0")
    assert backend == "scalar"


def test_backends_bit_identical():
    default_backend, default_digest = _run()
    if platform.machine().lower() in ("x86_64", "amd64") and default_backend != "avx":
        pytest.skip("CPU or OS without AVX: only the scalar path is available")
    scalar_backend, scalar_digest = _run("scalar")
    assert scalar_backend == "scalar"
    assert default_digest == scalar_digest
