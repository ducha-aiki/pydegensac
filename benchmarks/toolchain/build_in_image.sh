#!/usr/bin/env bash
# Build the pydegensac extension inside a manylinux image and drop the raw,
# *un-repaired* .so into /out.
#
# Un-repaired on purpose: auditwheel vendors the image's LAPACK into the wheel,
# which is a second variable. Copying the bare .so out lets every arm resolve
# LAPACK from the host, so the only thing left varying is the toolchain that
# compiled the code. openblas-devel is installed first for the same reason CI
# does it — CMake's FindLAPACK prefers OpenBLAS, so all arms link
# libopenblas.so.0 and pick up the host's copy at run time.
#
#   docker run --rm -v $PWD/..:/src:ro -v $OUT:/out <image> bash /src/... <tag>
set -euo pipefail

REF="${1:?usage: build_in_image.sh <git-ref>}"
PY=/opt/python/cp311-cp311/bin/python
# v_0.2.2 calls dgeqp3_ with no prototype in scope; GCC 14 (manylinux_2_28)
# rejects that outright, GCC 10 (manylinux2014) only warns. Demoting it back to
# a warning is what makes the two images comparable on the same old source —
# it is a diagnostic setting and changes no code generation. HEAD does not need
# it: the declaration was added in the lapwrap fix.
export CFLAGS="${EXTRA_CFLAGS:-} ${CFLAGS:-}"

echo "== image: $(cat /etc/*release | grep -m1 PRETTY_NAME || true)"
echo "== gcc: $(gcc --version | head -1)"

# CentOS 7 is EOL; manylinux2014 images ship vault-pinned repos, but be loud
# if that ever stops working rather than silently falling back to no LAPACK.
yum install -y openblas-devel lapack-devel blas-devel gcc-gfortran >/dev/null

# Clone rather than copy: /src carries a stale host CMakeCache under build/
# (which makes CMake refuse to configure) plus ~10 GB of benchmark data.
rm -rf /tmp/src
git config --global --add safe.directory '*'
git clone -q --no-local --no-checkout /src /tmp/src
cd /tmp/src
git checkout -q --detach "$REF"
git submodule update --init --recursive -q 2>/dev/null || true

rm -rf /tmp/pkg
# VERBOSE=1 makes the CMake-generated Makefiles echo every compile line, which
# is the evidence for "what flags did this toolchain actually use".
VERBOSE=1 "$PY" -m pip install --no-deps --target /tmp/pkg . > /out/build.log 2>&1 \
    || { tail -40 /out/build.log; exit 1; }

SO=$(find /tmp/pkg -name 'pydegensac*.so' | head -1)
[ -n "$SO" ] || { echo "no .so produced"; tail -40 /out/build.log; exit 1; }
cp "$SO" /out/
echo "== produced $(basename "$SO")"
{
    gcc --version | head -1
    cmake --version | head -1
    echo "dynsyms: $(nm -D --defined-only "$SO" | wc -l)"
    echo "-- .comment --"
    readelf -p .comment "$SO" || true
    echo "-- one exp_ranH.c compile line --"
    grep -m1 -o '/[^ ]*cc1\?[^ ]* .*exp_ranH\.c[^ ]*' /out/build.log || \
        grep -m1 'exp_ranH\.c' /out/build.log || true
} > /out/buildinfo.txt
cat /out/buildinfo.txt
