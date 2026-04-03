# Repository Guidelines

## Project Structure & Module Organization
`src/pydegensac/` contains the Python package, the `pybind11` binding layer in `bindings.cpp`, and the bundled C/C sources under `degensac/` and `matutls/`. Keep Python-facing helpers in `src/pydegensac/utils.py` and public exports in `src/pydegensac/__init__.py`. Put regression tests in `tests/`; the current suite is minimal and centered on import/build validation. Use `examples/` for runnable demos and notebooks, `docs/` for Sphinx docs, and treat `lib/pybind11/` as vendored third-party code unless a dependency update is intentional.

## Build, Test, and Development Commands
Install locally with `pip install .` for a normal build or `python setup.py build` while iterating on the extension. Run `pytest tests` to execute the repository test suite. Smoke-test the built module with `cd examples && python -utt simple-example.py`; this mirrors the wheel CI test command. For packaging work, GitHub Actions uses `python -m cibuildwheel --output-dir wheelhouse` from `.github/workflows/build_wheels.yml`.

## Coding Style & Naming Conventions
Follow the existing style instead of reformatting unrelated files. Python uses 4-space indentation, snake_case function names, and thin wrappers around the compiled extension. C/C++ code in `bindings.cpp` and the legacy sources uses concise, low-abstraction procedural code; preserve current naming and brace style when editing those files. No formatter or linter is configured in this repo, so keep diffs small and readable.

## Testing Guidelines
Add tests under `tests/` with `test_*.py` names and `pytest`-style assertions. Prefer small regression tests for Python API behavior, array shape validation, and import/build failures. When changing the compiled extension or build logic, also run the example smoke test so the extension is exercised end to end.

## Commit & Pull Request Guidelines
Recent history uses short, imperative commit subjects such as `fix the build` and `remove data copy`. Keep commit titles brief, specific, and focused on one change. Pull requests should describe the user-visible effect, note platform-specific build implications, and link related issues when relevant. Include example output or screenshots only when documentation or notebook content changes.

## Build Environment Notes
This project depends on CMake, BLAS/LAPACK, and a C++11-capable compiler. Cross-platform wheel builds are handled in GitHub Actions with `cibuildwheel`; update the workflow when changing binary dependencies or supported Python versions.
