#!/usr/bin/env bash
# Build a real, auditwheel-repaired wheel inside a manylinux image, the way CI
# does — and drop it in /out.
#
# The companion to build_in_image.sh, which deliberately copies the bare .so out
# so the compiler can be measured alone. This one keeps the other half: what
# `auditwheel repair` vendors into the wheel. A shipped wheel carries the
# image's own BLAS/LAPACK and uses it in preference to anything on the user's
# machine, so the image choice is two decisions, not one.
#
#   docker run --rm -v $PWD/..:/src:ro -v $OUT:/out <image> \
#       bash /src/benchmarks/toolchain/build_wheel_in_image.sh <git-ref>
set -euo pipefail

REF="${1:?usage: build_wheel_in_image.sh <git-ref>}"
PY=/opt/python/cp311-cp311/bin/python
export CFLAGS="${EXTRA_CFLAGS:-} ${CFLAGS:-}"

echo "== image: $(grep -m1 PRETTY_NAME /etc/os-release || true)"
echo "== gcc:   $(gcc --version | head -1)"

# Exactly CI's CIBW_BEFORE_ALL_LINUX, openblas-devel first so CMake's
# FindLAPACK prefers it over the reference build.
yum install -y openblas-devel lapack-devel blas-devel gcc-gfortran >/dev/null

rm -rf /tmp/src
git config --global --add safe.directory '*'
git clone -q --no-local --no-checkout /src /tmp/src
cd /tmp/src
git checkout -q --detach "$REF"

rm -rf /tmp/wh /tmp/rep
"$PY" -m pip wheel --no-deps -w /tmp/wh . > /out/build.log 2>&1 \
    || { tail -40 /out/build.log; exit 1; }
auditwheel repair -w /tmp/rep /tmp/wh/*.whl >> /out/build.log 2>&1 \
    || { tail -40 /out/build.log; exit 1; }
cp /tmp/rep/*.whl /out/

{
    echo "== $(basename /tmp/rep/*.whl)"
    echo "-- what auditwheel vendored --"
    "$PY" - <<'EOF'
import glob, zipfile
w = glob.glob("/tmp/rep/*.whl")[0]
for n in sorted(zipfile.ZipFile(w).namelist()):
    if ".libs/" in n and n.endswith(tuple(f".so.{i}" for i in range(10)) + (".so",)) or ".libs/" in n:
        print("  ", n.split("/")[-1])
EOF
    echo "-- openblas version in the image --"
    rpm -q openblas-devel lapack-devel 2>/dev/null || true
} > /out/wheelinfo.txt
cat /out/wheelinfo.txt
