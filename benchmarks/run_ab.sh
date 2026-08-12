#!/usr/bin/env bash
# Benchmark two pydegensac builds against each other and against the field.
#
#   ./run_ab.sh [BASE_REF]          # BASE_REF defaults to master
#   PYTHON=/path/to/python ./run_ab.sh
#
# Both builds are installed into their own directory with `pip --target`, and
# selected per run with PYTHONPATH — one interpreter, one set of cv2/poselib
# wheels, so the only thing that differs between the arms is pydegensac.
#
# $PYTHON must have numpy, opencv-python, poselib, h5py and matplotlib, plus a
# working build toolchain for the extension (a C++ compiler, CMake and LAPACK).
# See README.md if yours does not.
#
# Only the pydegensac import differs between the arms, so the base arm runs
# `--methods pydegensac` while the branch arm runs the whole roster.
# Results land in results/{f,h}_{base,branch}.jsonl; report.py merges them.
set -euo pipefail
cd "$(dirname "$0")"

PYTHON="${PYTHON:-python3}"
BASE_REF="${1:-master}"
BRANCH_REF="$(git rev-parse --abbrev-ref HEAD)"
REPO_ROOT="$(git rev-parse --show-toplevel)"
AB="$PWD/.ab"
WORKTREE="$AB/worktree-base"
mkdir -p "$AB"

# CMake and LAPACK often live in the interpreter's own prefix (conda), which
# is not on PATH unless the environment is activated — put it there for the
# build so `pip install .` finds both.
PY_PREFIX="$("$PYTHON" -c 'import sys; print(sys.prefix)')"
build() {  # name source_dir
    local name="$1" src="$2"
    rm -rf "$AB/pkg-$name"
    PATH="$PY_PREFIX/bin:$PATH" CMAKE_PREFIX_PATH="$PY_PREFIX" \
        "$PYTHON" -m pip install -q --no-deps --target "$AB/pkg-$name" "$src"
}

# The base arm builds from a detached worktree, so the working tree is never
# touched. Kept between runs; remove .ab to start clean.
if [ ! -d "$WORKTREE" ]; then
    git -C "$REPO_ROOT" worktree add --detach "$WORKTREE" "$BASE_REF"
fi
BASE_SHA="$(git -C "$WORKTREE" rev-parse --short HEAD)"
BRANCH_SHA="$(git -C "$REPO_ROOT" rev-parse --short HEAD)"

echo "== building pydegensac: base=$BASE_REF@$BASE_SHA branch=$BRANCH_REF@$BRANCH_SHA"
build base "$WORKTREE"
build branch "$REPO_ROOT"

run() {  # arm label extra-args...
    local arm="$1" label="$2"; shift 2
    PYTHONPATH="$AB/pkg-$arm" "$PYTHON" run.py "$@" --label "$label"
}

echo "== branch arm ($BRANCH_REF@$BRANCH_SHA): full roster"
run branch "branch:$BRANCH_REF@$BRANCH_SHA" f --out f_branch.jsonl
run branch "branch:$BRANCH_REF@$BRANCH_SHA" h --out h_branch.jsonl

echo "== base arm ($BASE_REF@$BASE_SHA): pydegensac only"
run base "base:$BASE_REF@$BASE_SHA" f --methods pydegensac --out f_base.jsonl
run base "base:$BASE_REF@$BASE_SHA" h --methods pydegensac --out h_base.jsonl

echo "== report"
PYTHONPATH="$AB/pkg-branch" "$PYTHON" report.py f \
    results/f_branch.jsonl results/f_base.jsonl --plot results/time_maa_f.png
PYTHONPATH="$AB/pkg-branch" "$PYTHON" report.py h \
    results/h_branch.jsonl results/h_base.jsonl --plot results/time_maa_h.png
