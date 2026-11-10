#!/usr/bin/env bash
# Build the current checkout and benchmark it against the full roster.
# Usage: PYTHON=/path/to/python ./run_current.sh [linux|m1|other-label]
set -euo pipefail
cd "$(dirname "$0")"

PYTHON="${PYTHON:-python3}"
PLATFORM="${1:-local}"
if [[ ! "$PLATFORM" =~ ^[a-zA-Z0-9_-]+$ ]]; then
    echo "Platform label must contain only letters, digits, underscores or hyphens." >&2
    exit 1
fi
REPO_ROOT="$(git rev-parse --show-toplevel)"
PKG="$PWD/.ab/pkg-current-$PLATFORM"
PY_PREFIX="$("$PYTHON" -c 'import sys; print(sys.prefix)')"
# Avoid stale dist-info from a previous release mislabeling the new binary.
rm -rf "$PKG"
PATH="$PY_PREFIX/bin:$PATH" CMAKE_PREFIX_PATH="$PY_PREFIX" \
    "$PYTHON" -m pip install --no-deps --target "$PKG" "$REPO_ROOT"
export PYTHONPATH="$PKG"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
VERSION="$("$PYTHON" -c 'from importlib.metadata import version; print(version("pydegensac"))')"
LABEL="$VERSION:$(git rev-parse --short HEAD)"

for problem in f h; do
    "$PYTHON" run.py "$problem" --label "$LABEL" --out "${problem}_current_${PLATFORM}.jsonl"
    "$PYTHON" report.py "$problem" "results/${problem}_current_${PLATFORM}.jsonl" \
        --full --plot "results/time_maa_${problem}_${PLATFORM}.png" \
        > "results/report_${problem}_${PLATFORM}.md"
done
