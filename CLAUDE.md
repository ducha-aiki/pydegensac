# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Python wrapper (via pybind11) around the original C implementation of LO-RANSAC / DEGENSAC (Chum et al.) for robust homography and fundamental matrix estimation. Published to PyPI as `pydegensac`.

## Build & test commands

Building requires CMake, LAPACK/BLAS, and a C++11 compiler. Use the platform default compiler — Clang on macOS. (The README's `CC=gcc-8` hint is a 2020 third-party note; current Clang builds fine, the 2026-08 profiling work was done under Clang, and the golden bit-exactness baselines were captured with it, so switching compilers risks perturbing FP codegen.)

```bash
pip install .                  # build + install (CMake is driven by setup.py)
python3 setup.py build         # build without installing, when iterating on the extension
pytest tests                   # run the test suite
pytest tests/test_ransac_smoke.py::test_public_api_seed_makes_results_deterministic  # single test
cd examples && python simple-example.py   # end-to-end smoke test (same command CI wheels use)
```

There is no linter or formatter configured; follow existing style and keep diffs small. When changing the compiled extension or build logic, also run the example smoke test so the extension is exercised end to end.

Wheels are built by `.github/workflows/build_wheels.yml` with `cibuildwheel` for Linux/macOS/Windows; update that workflow when changing binary dependencies or supported Python versions.

## Architecture

Three layers, from bottom up:

1. **Legacy C core** — `src/pydegensac/degensac/` (RANSAC/DEGENSAC algorithms, e.g. `exp_ranH.c`, `exp_ranF.c`, `ranH.c` (degeneracy-support subset), `DegUtils.c`, `Ftools.c`, `Htools.c`, `rtools.c`, `utools.c`, `hash.c`, `lapwrap.c`, `bsd_random.c` (lock-free BSD TYPE_3 RNG, bit-compatible with libc `random()` — see `docs/reports/2026-08-10-speedup-session.md`)) and `src/pydegensac/matutls/` (linear algebra utilities). This is decades-old procedural C linked against LAPACK; preserve its naming and style when editing. Compiled into static libs `pydegensac_support` and `matutls` by the top-level `CMakeLists.txt`. Golden fixed-seed regression tests in `tests/test_golden_regression.py` guarantee bit-equivalence; regenerate baselines only deliberately via `scripts/make_golden_data.py`.

2. **pybind11 binding** — `src/pydegensac/bindings.cpp` exposes the low-level functions `findHomography_` and `findFundamentalMatrix_` (note trailing underscore) that take all parameters positionally with error types as ints. `lib/pybind11/` is vendored third-party code — don't touch it unless a dependency update is intentional.

3. **Python public API** — `src/pydegensac/utils.py` defines `findHomography` / `findFundamentalMatrix`, thin wrappers that validate inputs, map error-type strings ("sampson", etc.) to the ints the C layer expects, and handle LAF consistency options. Exported via `src/pydegensac/__init__.py`.

Build plumbing: `setup.py` defines a custom `CMakeBuild` command that invokes CMake (including per-arch macOS handling for cibuildwheel); `CMakeLists.txt` ties the three layers together.

### Non-obvious behaviors in the Python layer

- Input keypoints may be numpy arrays of shape Nx2 (x, y) or Nx6 (x, y + flattened 2x2 affine frame), or a list of `cv2.KeyPoint` (converted via `convert_cv2_kpts_to_xyA`). The Nx6 form enables the LAF consistency check (`laf_consistensy_coef`); with Nx2 input it is silently disabled with a warning.
- `findHomography` post-processes the matrix returned from C++ with `np.linalg.inv(H.T)` — the C layer works with the transposed inverse convention.
- An all-zeros model from the C layer means "no good model found"; the wrappers then return an all-`False` mask.
- Both estimators accept `seed` (default -1 = nondeterministic); a fixed seed makes results reproducible, covered by `tests/test_ransac_smoke.py`.

## Tests

Tests live in `tests/` as `test_*.py` with pytest-style assertions. The suite is small: import validation plus smoke tests that exercise both the public API and the underscore low-level bindings with synthetic point correspondences. Prefer small regression tests for Python API behavior and array-shape validation.

## Commits & PRs

Short, imperative commit subjects (e.g. "fix the build"). PRs should describe the user-visible effect and note platform-specific build implications.
