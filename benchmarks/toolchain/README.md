# Which toolchain built this wheel, and what did it cost?

Published Linux wheels used to run 1.10x (F) / 1.30x (H) slower than the same
source built locally. This harness is what settled why: it rebuilds one git ref
inside each manylinux image, pulls the bare `.so` out, and times them all
against a local build on one core with a fixed seed.

Answer, and the reason CI builds on `manylinux_2_28`: it was the compiler, and
almost all of it was one function, `pinvJ`. Full write-up in
`docs/reports/2026-08-11-ransac-benchmark.md` ("The wheel gap is `pinvJ`").

## Why it is built this way

- **The `.so` is copied out un-repaired**, not the wheel. `auditwheel` vendors
  the image's LAPACK into the wheel, so a wheel-vs-wheel comparison moves two
  variables at once. A bare `.so` resolves `libopenblas.so.0` from the host, so
  every arm shares one LAPACK and only the compiler differs.
- **`openblas-devel` is installed in the image first**, as CI does, so CMake's
  `FindLAPACK` picks OpenBLAS in every image and the soname matches the host's.
- **Fixed seed.** Every arm does bit-identical work; `summarize.py` prints the
  inlier checksum and says whether the arms agree. If they don't, the timings
  are meaningless.
- **Interleaved rounds, pinned to a P-core.** On a hybrid CPU, drifting onto an
  E-core costs ~25% — enough to invent a difference that isn't there.
- **A dedicated driver, not `run.py`.** The benchmark's reprojection metric is
  85% of its process time and buries the estimator.

## Running it

```bash
cd benchmarks
python toolchain/dump_pairs.py          # once; needs the datasets (setup_data.py)

# build one ref in one image; the .so lands in .ab/iso/out-<name>/
docker run --rm \
    -v "$PWD/..":/src:ro -v "$PWD/.ab/iso/out-ml2014":/out \
    quay.io/pypa/manylinux2014_x86_64 \
    bash /src/benchmarks/toolchain/build_in_image.sh HEAD

# assemble an arm: a local install with the container's .so swapped in
cp -r .ab/pkg-local-head .ab/pkg-ml2014-head
cp .ab/iso/out-ml2014/pydegensac*.so .ab/pkg-ml2014-head/pydegensac/

ARMS="pkg-local-head pkg-ml2014-head" toolchain/run_arms.sh h 7 > h.jsonl
python toolchain/summarize.py h.jsonl
```

`build_wheel_in_image.sh` is the companion for the other half of the question:
it runs CI's `CIBW_BEFORE_ALL_LINUX` and `auditwheel repair`, so you get the
wheel *as shipped*, vendored BLAS and all. Use `build_in_image.sh` to compare
compilers (one LAPACK, held constant) and this one to compare images as users
receive them.

```bash
docker run --rm -v "$PWD/..":/src:ro -v "$PWD/.ab/iso/whl-ml228":/out \
    quay.io/pypa/manylinux_2_28_x86_64 \
    bash /src/benchmarks/toolchain/build_wheel_in_image.sh HEAD
pip install --no-deps --target .ab/pkg-whl-ml228 .ab/iso/whl-ml228/*.whl
```

`EXTRA_CFLAGS` is passed through to the container build — that is how the
visibility hypothesis was tested, and how pre-`lapwrap` refs are built at all
(`v_0.2.2` needs `-Wno-error=implicit-function-declaration`, since gcc 14
rejects its undeclared `dgeqp3_`).

To attribute a difference to a function:

```bash
LD_LIBRARY_PATH=$HOME/miniconda3/envs/py313/lib \
PYTHONPATH=$PWD/.ab/pkg-ml2014-head \
    perf record -F 2000 -g -o /tmp/p.data -- \
    taskset -c 2 .ab/env_np1/bin/python toolchain/iso_bench.py h --reps 1
perf report -i /tmp/p.data --no-children --stdio -F sample,symbol -g none
```

WSL2 exposes no hardware PMU, so this is timer sampling only — good enough to
attribute time to a symbol, useless for instruction counts or IPC.
