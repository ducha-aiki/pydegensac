#!/usr/bin/env bash
# Compare published PyPI releases against local builds of the same code, to
# answer: did anything about the Linux build get slower along the way?
#
#   PYTHON=/path/to/python ./run_releases.sh
#
# A PyPI wheel differs from a local build in two ways at once — the source, and
# the build environment (manylinux gcc 10.2 + reference LAPACK from
# `yum lapack-devel`, versus whatever your box has). So this runs the 0.2.2
# *tag* built locally alongside the 0.2.2 *wheel*: the wheel-vs-tag gap is the
# build environment, the tag-vs-master gap is the source.
#
# $PYTHON must be numpy<2: releases 0.1.2 and 0.2 were built against an old
# pybind11 and silently return every point as an inlier under numpy 2.x
# (0.2 is yanked for this; 0.1.2 is not). Holding numpy fixed across all arms
# also keeps the release axis clean.
set -euo pipefail
cd "$(dirname "$0")"

PYTHON="${PYTHON:-$PWD/.ab/env_np1/bin/python}"
AB="$PWD/.ab"
REPO_ROOT="$(git rev-parse --show-toplevel)"
PY_PREFIX="$("$PYTHON" -c 'import sys; print(sys.prefix)')"

"$PYTHON" - <<'EOF'
import numpy
assert numpy.__version__ < "2", (
    f"numpy {numpy.__version__}: releases 0.1.2/0.2 return garbage under "
    "numpy 2.x. Use a numpy<2 interpreter (see the header of this script).")
EOF

#: wheels straight from PyPI
PYPI_VERSIONS="${PYPI_VERSIONS:-0.1.2 0.2.1 0.2.2}"
#: git refs built here, in this machine's toolchain
LOCAL_REFS="${LOCAL_REFS:-v_0.2.2 master}"

build_local() {  # name ref
    local name="$1"
    local ref="$2"
    local wt="$AB/wt-$name"
    [ -d "$wt" ] || git -C "$REPO_ROOT" worktree add --detach "$wt" "$ref"
    rm -rf "$AB/pkg-$name"
    PATH="$PY_PREFIX/bin:$PATH" CMAKE_PREFIX_PATH="$PY_PREFIX" \
        "$PYTHON" -m pip install -q --no-deps --target "$AB/pkg-$name" "$wt"
}

arm() {  # pkgdir label
    local pkg="$1"
    local label="$2"
    local tag="${label%%:*}"
    echo "== $label"
    PYTHONPATH="$pkg" "$PYTHON" run.py f --methods pydegensac \
        --label "$label" --out "rel_f_$tag.jsonl"
    PYTHONPATH="$pkg" "$PYTHON" run.py h --methods pydegensac \
        --label "$label" --out "rel_h_$tag.jsonl"
}

for v in $PYPI_VERSIONS; do
    rm -rf "$AB/pkg-pypi-$v"
    "$PYTHON" -m pip install -q --no-deps --target "$AB/pkg-pypi-$v" \
        "pydegensac==$v" 2>/dev/null
done
for ref in $LOCAL_REFS; do
    build_local "local-${ref#v_}" "$ref"
done
rm -rf "$AB/pkg-local-branch"
PATH="$PY_PREFIX/bin:$PATH" CMAKE_PREFIX_PATH="$PY_PREFIX" \
    "$PYTHON" -m pip install -q --no-deps --target "$AB/pkg-local-branch" \
    "$REPO_ROOT"

for v in $PYPI_VERSIONS; do
    arm "$AB/pkg-pypi-$v" "pypi-$v:wheel"
done
for ref in $LOCAL_REFS; do
    n="local-${ref#v_}"
    arm "$AB/pkg-$n" "$n:$(git -C "$AB/wt-$n" rev-parse --short HEAD)"
done
arm "$AB/pkg-local-branch" "local-branch:$(git -C "$REPO_ROOT" rev-parse --short HEAD)"

echo "== report"
"$PYTHON" report.py f results/rel_f_*.jsonl --full \
    --plot results/releases_f.png > results/report_releases_f.md
"$PYTHON" report.py h results/rel_h_*.jsonl --full \
    --plot results/releases_h.png > results/report_releases_h.md
echo "wrote results/report_releases_{f,h}.md and results/releases_{f,h}.png"
