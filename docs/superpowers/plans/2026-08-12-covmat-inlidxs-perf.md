# `cov_mat` / `inlidxs` Performance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the two structural hot spots the post-LAPACK-fix profile exposes — `cov_mat`'s 45 strided passes (60% of H) and `inlidxs`'s per-point cross-TU division (46% of F) — without changing the algorithm.

**Architecture:** Four code changes, ordered so the two bit-exact ones land first and the two numerics-changing ones land separately behind a statistical gate. A new committed harness (`scripts/stat_ab.py`) provides that gate, replacing the ad hoc protocol used for `4430bc7`.

**Tech Stack:** C99 (the legacy core, built by CMake into `pydegensac_support`), pybind11 binding, pytest for the golden gate, numpy/scipy for the statistical harness.

## Global Constraints

- **No algorithmic change.** Same models, same iteration sequence, same accept/reject decisions. Only floating-point rounding may differ.
- **No SIMD intrinsics, no architecture-specific code.** Wheels are portable; changes must be plain C99.
- **`cov_mat` must stay correct for arbitrary `siz`**, not just 9 — `exp_ranH.c:750` passes `nullsize`.
- **Build with the platform default compiler (Clang on macOS)**, per CLAUDE.md. Do not switch compilers; it perturbs FP codegen and the golden baselines.
- After any change to the compiled extension: `pip install .`, then `pytest tests`, then `cd examples && python3 simple-example.py`.
- The golden gate runs in exact mode locally (`PYDEGENSAC_GOLDEN_EXACT=1`, the default).
- Statistical gate threshold: the `4430bc7` run's minimum KS p-value was **0.108**. A comparable minimum passes; a small p is stop-and-investigate, not a threshold to tune.
- Measurements are **M1 only** in this session. Label them as such; the maintainer re-runs on Linux separately.

---

### Task 1: Statistical A/B harness

**Files:**
- Create: `scripts/stat_ab.py`

**Interfaces:**
- Consumes: nothing from earlier tasks. Reads the golden `.npz` fixtures in `tests/data/` (keys `pts1`, `pts2`, `px_th`, `conf`, `max_iters`, plus `H_gt` for H pairs and `F_gt`/`K1`/`K2`/`R1`/`R2`/`T1`/`T2` for F pairs).
- Produces: CLI `python3 scripts/stat_ab.py --build-a DIR --build-b DIR [--runs N]`, where each `DIR` is a `pip --target` directory containing `pydegensac`. Prints one row per (pair, metric) with the KS statistic and p-value, then a final `min p = X` line. Exit code 1 if any p < 0.01.

- [ ] **Step 1: Write the harness**

Create `scripts/stat_ab.py`:

```python
#!/usr/bin/env python3
"""Statistical A/B for deliberate numerics changes.

Runs both builds over the golden pairs at many seeds and compares the
resulting distributions. Used to gate changes that alter floating-point
results without altering the algorithm -- the protocol first used for the
per-iteration reseed removal (4430bc7), kept here so it is reproducible.

    python3 scripts/stat_ab.py --build-a .ab/pkg-base --build-b .ab/pkg-branch

Each build directory is a `pip install --target` tree containing pydegensac.
The two builds run in separate subprocesses, because pydegensac can only be
imported once per process.
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
DATA_DIR = REPO / "tests" / "data"
GT_INLIER_PX = 3.0

# Worker: runs inside a subprocess with one build on sys.path, writes JSON to
# stdout. Kept as a string so the harness stays a single file.
WORKER = r"""
import json, sys
import numpy as np
sys.path.insert(0, sys.argv[1])
import pydegensac

path, runs = sys.argv[2], int(sys.argv[3])
d = np.load(path)
pts1, pts2 = d["pts1"], d["pts2"]
is_h = "H_gt" in d

def gt_err(model, mask):
    p1 = np.concatenate([pts1[mask], np.ones((mask.sum(), 1))], axis=1)
    if is_h:
        proj = (d["H_gt"] @ p1.T).T
        return np.linalg.norm(proj[:, :2] / proj[:, 2:3] - pts2[mask], axis=1)
    p2 = np.concatenate([pts2[mask], np.ones((mask.sum(), 1))], axis=1)
    Fx1 = (d["F_gt"] @ p1.T).T
    num = np.abs(np.sum(p2 * Fx1, axis=1))
    return num / np.sqrt(Fx1[:, 0] ** 2 + Fx1[:, 1] ** 2)

# GT-consistent correspondences, for the precision metric
all_idx = np.ones(len(pts1), bool)
gt_inlier = gt_err(None, all_idx) <= 3.0

out = {"inliers": [], "gt_err": [], "gt_prec": []}
for seed in range(1, runs + 1):
    if is_h:
        M, mask = pydegensac.findHomography(
            pts1, pts2, px_th=float(d["px_th"]), conf=float(d["conf"]),
            max_iters=int(d["max_iters"]), seed=seed)
    else:
        M, mask = pydegensac.findFundamentalMatrix(
            pts1, pts2, px_th=float(d["px_th"]), conf=float(d["conf"]),
            max_iters=int(d["max_iters"]), seed=seed)
    mask = np.asarray(mask, bool)
    n = int(mask.sum())
    out["inliers"].append(n)
    out["gt_err"].append(float(np.median(gt_err(M, mask))) if n else float("nan"))
    out["gt_prec"].append(float(gt_inlier[mask].mean()) if n else float("nan"))
print(json.dumps(out))
"""


def run_build(build_dir, npz_path, runs):
    proc = subprocess.run(
        [sys.executable, "-c", WORKER, str(build_dir), str(npz_path), str(runs)],
        capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"worker failed for {npz_path.name}:\n{proc.stderr}")
    return json.loads(proc.stdout)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--build-a", required=True, help="reference build directory")
    ap.add_argument("--build-b", required=True, help="candidate build directory")
    ap.add_argument("--runs", type=int, default=1000, help="seeds per pair")
    args = ap.parse_args()

    pairs = sorted(DATA_DIR.glob("golden_*.npz"))
    if not pairs:
        raise SystemExit(f"no golden fixtures in {DATA_DIR}")

    print(f"{'pair':52s} {'metric':9s} {'KS':>8} {'p':>10}")
    worst = (1.0, None)
    for npz in pairs:
        a = run_build(args.build_a, npz, args.runs)
        b = run_build(args.build_b, npz, args.runs)
        for metric in ("inliers", "gt_err", "gt_prec"):
            xa = np.asarray(a[metric], float)
            xb = np.asarray(b[metric], float)
            xa, xb = xa[np.isfinite(xa)], xb[np.isfinite(xb)]
            ks, p = stats.ks_2samp(xa, xb)
            print(f"{npz.stem[:52]:52s} {metric:9s} {ks:8.4f} {p:10.4g}")
            if p < worst[0]:
                worst = (p, f"{npz.stem} / {metric}")
    print(f"\nmin p = {worst[0]:.4g}  ({worst[1]})   over "
          f"{len(pairs) * 3} comparisons, {args.runs} seeds each")
    if worst[0] < 0.01:
        print("FAIL: distributions differ")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 2: Build two identical trees to calibrate the harness**

```bash
cd /Users/oldufo/dev/pydegensac
PY=benchmarks/.ab/env/bin/python
PATH="$PWD/benchmarks/.ab/env/bin:$PATH" CMAKE_PREFIX_PATH="$PWD/benchmarks/.ab/env" \
    $PY -m pip install -q --no-deps --target /tmp/sab-a .
cp -r /tmp/sab-a /tmp/sab-b
benchmarks/.ab/env/bin/pip install -q scipy
```

- [ ] **Step 3: Run it on identical builds — must pass**

Run: `benchmarks/.ab/env/bin/python scripts/stat_ab.py --build-a /tmp/sab-a --build-b /tmp/sab-b --runs 200`
Expected: every p large, final line `min p = ...` well above 0.01, exit 0. (Identical builds at identical seeds produce identical output, so KS p should be 1.0 throughout.)

- [ ] **Step 4: Verify it can detect a real difference**

The `master` build differs from the branch (RNG overhaul), so the harness must flag it.

```bash
benchmarks/.ab/env/bin/python scripts/stat_ab.py \
    --build-a benchmarks/.ab/pkg-base --build-b /tmp/sab-b --runs 200
```

Expected: differences visible. This is a *sensitivity check on the harness*, not a claim about `master` — the reseed removal was itself validated as distributionally identical, so p-values may well be large; what matters is that the run completes and produces per-pair statistics from two genuinely different binaries.

- [ ] **Step 5: Commit**

```bash
git add scripts/stat_ab.py
git commit -m "add the statistical A/B harness as a script

The 1000-seed KS protocol used to gate the reseed removal was ad hoc and not
kept. This is the second deliberate numerics change, so it becomes a committed
tool: two pip --target builds in, per-pair KS tests on inlier count, GT error
and GT precision out."
```

---

### Task 2: Inline `truncQuad` and make the compress store branchless (bit-exact)

**Files:**
- Modify: `src/pydegensac/degensac/rtools.h:74` (declaration → `static inline` definition)
- Modify: `src/pydegensac/degensac/rtools.c:163-174` (`inlidxs`), `rtools.c:231-239` (remove the out-of-line definition)

**Interfaces:**
- Consumes: `Score` struct from `rtools.h` (fields `I`, `J`, `Is`, `Il`).
- Produces: `truncQuad(double epsilon, double thr)` as a `static inline` function in `rtools.h`, same signature and same value as before. Callers in `exp_ranF.c:487,637` and `exp_ranH.c:436,559` keep working unchanged.

- [ ] **Step 1: Establish the baseline the test compares against**

The test for this task is the golden gate in exact mode — this change must not move a single bit.

Run: `python3 -m pytest tests -q`
Expected: `33 passed`. If it is not green before you start, stop and fix that first.

- [ ] **Step 2: Move `truncQuad` into the header as `static inline`**

In `src/pydegensac/degensac/rtools.h`, replace the declaration on line 74:

```c
double truncQuad(double epsilon, double thr);
```

with the definition:

```c
/* Truncated-quadratic gain. Defined here rather than in rtools.c so it
   inlines: it is called once per correspondence in inlidxs's inner loop, and
   a cross-TU call there was measurable (46% of F runtime sat in that loop). */
static inline double truncQuad(double epsilon, double thr) {
  if (thr == 0) {
      return 0;
    }
  if ( epsilon >= thr*9/4 ) {
      return 0;
    }
  return 1 - (epsilon/(thr*9/4));
}
```

Then delete the out-of-line definition at `rtools.c:231-239`.

- [ ] **Step 3: Make the compress store branchless in `inlidxs`**

In `src/pydegensac/degensac/rtools.c`, replace the body of `inlidxs`:

```c
Score inlidxs (const double * err, int len, double th, int * inl) {
  unsigned i;
  Score s = {0,0,0,0};
  for (i = 0; i < len; ++i) {
      s.J += truncQuad(err[i], th);
      if (err[i] <= th) {
          inl[s.I] = i;
          ++(s.I);
        }
    }
  return s;
}
```

with:

```c
Score inlidxs (const double * err, int len, double th, int * inl) {
  int i;
  Score s = {0,0,0,0};
  /* Branchless compress store: the write always happens, the index only
     advances for inliers. `inl` is allocated at `len` entries by every caller,
     so the speculative write at s.I is always in bounds. */
  for (i = 0; i < len; ++i) {
      s.J += truncQuad(err[i], th);
      inl[s.I] = i;
      s.I += (err[i] <= th);
    }
  return s;
}
```

- [ ] **Step 4: Verify bit-exactness**

```bash
pip install -q .
python3 -m pytest tests -q
```

Expected: `33 passed`. Any failure means the change was not bit-exact — revert and investigate before continuing; do not regenerate baselines for this task.

- [ ] **Step 5: Measure**

```bash
cd /private/tmp/claude-501/-Users-oldufo-dev-pydegensac/01ca7b38-2a29-4465-b461-8eda4565f838/scratchpad
python3 profile_driver.py f 20 && python3 profile_driver.py h 20
```

Record the calls/second from each line of output, before and after, in the commit message.

- [ ] **Step 6: Commit**

```bash
git add src/pydegensac/degensac/rtools.h src/pydegensac/degensac/rtools.c
git commit -m "inline truncQuad and drop the branch in inlidxs

Both bit-exact: the golden gate stays green in exact mode. truncQuad was a
cross-TU call per correspondence in the hottest loop in F (46% self time), and
the inlier test was a data-dependent branch."
```

---

### Task 3: Hoist the reciprocal out of the `inlidxs` loop (numerics change)

**Files:**
- Modify: `src/pydegensac/degensac/rtools.c` (`inlidxs`)
- Modify: `src/pydegensac/degensac/rtools.h` (`truncQuad`, add a scaled variant)

**Interfaces:**
- Consumes: `truncQuad` from Task 2.
- Produces: `truncQuadInv(double epsilon, double inv)` as `static inline` in `rtools.h`, where `inv == 1.0/(thr*9/4)`; returns `fmax(0.0, 1.0 - epsilon*inv)`. `truncQuad` itself is unchanged, so the four call sites in `exp_ranF.c`/`exp_ranH.c` are untouched.

- [ ] **Step 1: Add the scaled variant to `rtools.h`**

Add below `truncQuad`, and add `#include <math.h>` to `rtools.h` if it is not already included:

```c
/* truncQuad with the reciprocal pre-computed by the caller: inv = 1/(thr*9/4).
   Equivalent in exact arithmetic -- truncQuad returns 0 exactly where
   1 - epsilon*inv <= 0 -- but x*(1/y) does not round like x/y, so this is a
   deliberate numerics change (see docs/superpowers/specs/
   2026-08-12-covmat-inlidxs-perf-design.md). Hoisting the division out of the
   per-correspondence loop is the point. */
static inline double truncQuadInv(double epsilon, double inv) {
  return fmax(0.0, 1.0 - epsilon*inv);
}
```

- [ ] **Step 2: Use it in `inlidxs`**

Careful with `th == 0`: the original `truncQuad` returns 0 for it, but
`inv = 1/0 = inf` would make `truncQuadInv` return 0 only for `epsilon > 0` and
`fmax(0, NaN)` for `epsilon == 0`. Keep the `th == 0` case out of the loop with
a loop-invariant flag the compiler can hoist, rather than folding it into `inv`:

```c
Score inlidxs (const double * err, int len, double th, int * inl) {
  int i;
  Score s = {0,0,0,0};
  const double inv = 1.0/(th*9/4);
  const int score_j = (th != 0);
  for (i = 0; i < len; ++i) {
      if (score_j) s.J += truncQuadInv(err[i], inv);
      inl[s.I] = i;
      s.I += (err[i] <= th);
    }
  return s;
}
```

`score_j` is loop-invariant, so the compiler hoists the branch out of the loop; correctness for `th == 0` is preserved.

- [ ] **Step 3: Build and confirm the gate now fails**

```bash
pip install -q .
python3 -m pytest tests -q
```

Expected: **failures** in `test_*_bit_equivalent`. That is the expected consequence of a numerics change — it confirms the change is live. Do not regenerate baselines yet; that happens once, in Task 6.

- [ ] **Step 4: Run the statistical gate**

```bash
PATH="$PWD/benchmarks/.ab/env/bin:$PATH" CMAKE_PREFIX_PATH="$PWD/benchmarks/.ab/env" \
    benchmarks/.ab/env/bin/python -m pip install -q --no-deps --target /tmp/sab-recip .
benchmarks/.ab/env/bin/python scripts/stat_ab.py \
    --build-a /tmp/sab-a --build-b /tmp/sab-recip --runs 1000
```

Expected: exit 0, `min p` comparable to the 0.108 precedent. If any p is below 0.01, stop — the change is not distribution-preserving and needs investigation, not a looser threshold.

- [ ] **Step 5: Measure and commit**

```bash
cd /private/tmp/claude-501/-Users-oldufo-dev-pydegensac/01ca7b38-2a29-4465-b461-8eda4565f838/scratchpad
python3 profile_driver.py f 20
cd /Users/oldufo/dev/pydegensac
git add src/pydegensac/degensac/rtools.h src/pydegensac/degensac/rtools.c
git commit -m "hoist the per-correspondence division out of inlidxs

<calls/sec before> -> <calls/sec after> on F. Numerics change: x*(1/y) does
not round like x/y. Statistical gate over 1000 seeds x 10 golden pairs x 3
metrics: min p = <value>."
```

---

### Task 4: `cov_mat` in one pass (numerics change)

**Files:**
- Modify: `src/pydegensac/degensac/utools.c:170-184`
- Create: `tests/c/test_cov_mat.c` (correctness harness, compiled ad hoc; binaries land in the gitignored `tests/bin/`)

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `void cov_mat(double *Cv, const double *Z, int len, int siz)` — unchanged signature, unchanged semantics (`Cv = Zᵀ Z`, both triangles filled). All five call sites (`Ftools.c:367,422`, `Htools.c:125`, `exp_ranH.c:750,756`) are untouched.

- [ ] **Step 1: Write the failing test**

Create `tests/c/test_cov_mat.c`:

```c
/* Checks cov_mat against a naive reference over many shapes.
   Build and run:
     mkdir -p tests/bin && cc -O2 -I src/pydegensac/degensac \
       tests/c/test_cov_mat.c src/pydegensac/degensac/utools.c \
       -o tests/bin/test_cov_mat && ./tests/bin/test_cov_mat            */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void cov_mat(double *Cv, const double * Z, int len, int siz);

static void reference(double *Cv, const double *Z, int len, int siz) {
    int i, j, k;
    for (i = 0; i < siz; i++)
        for (j = 0; j <= i; j++) {
            double val = 0;
            for (k = 0; k < len*siz; k += siz) val += Z[k+i] * Z[k+j];
            Cv[siz*i + j] = val;
            Cv[i + siz*j] = val;
        }
}

int main(void) {
    const int sizes[] = {2, 3, 7, 8, 9, 10};
    const int lens[] = {1, 2, 8, 10, 33, 64, 257};
    int si, li, t, i;
    double worst = 0.0;
    srand(1);
    for (si = 0; si < 6; si++)
        for (li = 0; li < 7; li++)
            for (t = 0; t < 5; t++) {
                int siz = sizes[si], len = lens[li];
                double *Z = malloc(sizeof(double)*len*siz);
                double *A = malloc(sizeof(double)*siz*siz);
                double *B = malloc(sizeof(double)*siz*siz);
                for (i = 0; i < len*siz; i++)
                    Z[i] = (double)rand()/RAND_MAX*2.0 - 1.0;
                reference(A, Z, len, siz);
                cov_mat(B, Z, len, siz);
                for (i = 0; i < siz*siz; i++) {
                    double denom = fabs(A[i]) > 1.0 ? fabs(A[i]) : 1.0;
                    double rel = fabs(A[i] - B[i]) / denom;
                    if (rel > worst) worst = rel;
                }
                free(Z); free(A); free(B);
            }
    printf("worst relative difference vs reference: %.3g\n", worst);
    if (worst > 1e-12) { printf("FAIL\n"); return 1; }
    printf("PASS\n");
    return 0;
}
```

- [ ] **Step 2: Run it against the current implementation**

```bash
mkdir -p tests/bin
cc -O2 -I src/pydegensac/degensac tests/c/test_cov_mat.c \
   src/pydegensac/degensac/utools.c -o tests/bin/test_cov_mat && ./tests/bin/test_cov_mat
```

Expected: `PASS` with worst difference `0` — the test currently compares the implementation against a copy of itself. This proves the harness runs; it starts constraining behaviour the moment Step 3 changes the implementation.

- [ ] **Step 3: Rewrite `cov_mat` as a single pass**

Replace `utools.c:170-184` with:

```c
/* Cv = Z^T Z for Z of shape len x siz, row-major; both triangles filled.

   One pass over the points, accumulating every unique entry at once. The
   previous form ran one pass per entry -- 45 of them at siz=9 -- and each was
   a single serial FP accumulation chain, so it was latency-bound rather than
   throughput-bound. The hot caller (exp_ranH's per-iteration MCE solve, len=10)
   ran this on every RANSAC iteration.

   Summation order over points is unchanged; the accumulators are independent,
   which is what lets the FMA units fill. Rounding may differ from the old form
   in the last bits -- deliberate, see docs/superpowers/specs/
   2026-08-12-covmat-inlidxs-perf-design.md. */
void cov_mat(double *Cv, const double * Z, int len, int siz)
{
   int i, j, k;

   for (i = 0; i < siz*siz; i++) Cv[i] = 0.0;

   for (k = 0; k < len*siz; k += siz)
      for (i = 0; i < siz; i++) {
         const double zi = Z[k+i];
         for (j = 0; j <= i; j++)
            Cv[siz*i + j] += zi * Z[k+j];
      }

   for (i = 0; i < siz; i++)
      for (j = 0; j < i; j++)
         Cv[i + siz*j] = Cv[siz*i + j];
}
```

- [ ] **Step 4: Run the C test**

```bash
cc -O2 -I src/pydegensac/degensac tests/c/test_cov_mat.c \
   src/pydegensac/degensac/utools.c -o tests/bin/test_cov_mat && ./tests/bin/test_cov_mat
```

Expected: `PASS`, worst relative difference small but now possibly non-zero (different accumulation order). If it exceeds `1e-12`, the implementation is wrong, not merely differently rounded.

- [ ] **Step 5: Build, confirm the golden gate moves, run the statistical gate**

```bash
pip install -q .
python3 -m pytest tests -q      # expect bit-equivalence failures: numerics changed
PATH="$PWD/benchmarks/.ab/env/bin:$PATH" CMAKE_PREFIX_PATH="$PWD/benchmarks/.ab/env" \
    benchmarks/.ab/env/bin/python -m pip install -q --no-deps --target /tmp/sab-cov .
benchmarks/.ab/env/bin/python scripts/stat_ab.py \
    --build-a /tmp/sab-recip --build-b /tmp/sab-cov --runs 1000
```

Expected: statistical gate exit 0, `min p` comparable to 0.108.

- [ ] **Step 6: Measure and commit**

```bash
cd /private/tmp/claude-501/-Users-oldufo-dev-pydegensac/01ca7b38-2a29-4465-b461-8eda4565f838/scratchpad
python3 profile_driver.py h 20 && python3 profile_driver.py f 20
cd /Users/oldufo/dev/pydegensac
git add src/pydegensac/degensac/utools.c tests/c/test_cov_mat.c
git commit -m "compute cov_mat in one pass instead of one per entry

<calls/sec before> -> <calls/sec after> on H. 45 strided passes, each a serial
FP chain, become one pass with independent accumulators. Numerics change;
statistical gate min p = <value>. tests/c/test_cov_mat.c checks the result
against a naive reference across siz 2..10 and len 1..257."
```

---

### Task 5: Decide on `dsyrk` for the large-`len` call sites

**Files:**
- Modify: `src/pydegensac/degensac/utools.c` (only if the measurement says yes)

**Interfaces:**
- Consumes: `cov_mat` from Task 4.
- Produces: either nothing (measurement recorded, no code change) or a length-thresholded `dsyrk_` path inside `cov_mat` with an unchanged signature.

- [ ] **Step 1: Measure whether BLAS beats the single pass at large `len`**

Write `/tmp/bench_covmat.c` comparing the Task 4 implementation against `dsyrk_` for `siz=9`, `len` in {8, 10, 64, 256, 1024, 4096}, timing many repetitions of each:

```c
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
void cov_mat(double *Cv, const double * Z, int len, int siz);
extern void dsyrk_(char*, char*, ptrdiff_t*, ptrdiff_t*, double*, double*,
                   ptrdiff_t*, double*, double*, ptrdiff_t*);

int main(void) {
    const int lens[] = {8, 10, 64, 256, 1024, 4096};
    for (int li = 0; li < 6; li++) {
        int len = lens[li], siz = 9, reps = 200000000/(len*45) + 10;
        double *Z = malloc(sizeof(double)*len*siz), C[81];
        for (int i = 0; i < len*siz; i++) Z[i] = (double)rand()/RAND_MAX;
        clock_t t0 = clock();
        for (int r = 0; r < reps; r++) cov_mat(C, Z, len, siz);
        double own = (double)(clock()-t0)/CLOCKS_PER_SEC/reps*1e9;
        ptrdiff_t n = siz, k = len, lda = siz, ldc = siz;
        double one = 1.0, zero = 0.0;
        t0 = clock();
        for (int r = 0; r < reps; r++)
            dsyrk_("L", "N", &n, &k, &one, Z, &lda, &zero, C, &ldc);
        double blas = (double)(clock()-t0)/CLOCKS_PER_SEC/reps*1e9;
        printf("len %5d: single-pass %8.1f ns   dsyrk %8.1f ns\n", len, own, blas);
        free(Z);
    }
    return 0;
}
```

```bash
cc -O2 -I src/pydegensac/degensac /tmp/bench_covmat.c \
   src/pydegensac/degensac/utools.c -framework Accelerate -o tests/bin/bench_covmat
./tests/bin/bench_covmat
```

- [ ] **Step 2: Decide**

Adopt `dsyrk` only if it is faster by a clear margin (>20%) at a `len` that the estimators actually reach, and only above a threshold where it wins. `Ftools.c:367,422` and `Htools.c:125` pass `len` = inlier count, so those are the sites that could benefit; `exp_ranH.c:750,756` pass 8 and 10 and must keep the single-pass path.

If adopting, note that Fortran reads the row-major `len x siz` array `Z` as a column-major `siz x len` matrix, so the call is `dsyrk("L", "N", siz, len, 1.0, Z, siz, 0.0, Cv, siz)`, and the opposite triangle still needs filling.

- [ ] **Step 3: Record the outcome**

Either commit the thresholded implementation with the measurement table in the message, or — if BLAS does not win — record the numbers in `docs/reports/2026-08-12-benchmark-handoff.md` under "Smaller items" so nobody re-tests it:

```bash
git commit -am "record the dsyrk measurement: <verdict>"
```

---

### Task 6: Regenerate baselines, re-run the benchmark, update the write-ups

**Files:**
- Modify: `tests/data/golden_*.npz` (regenerated)
- Modify: `docs/reports/2026-08-11-ransac-benchmark.md` (new subsection under the M1 section)
- Modify: `docs/reports/2026-08-12-benchmark-handoff.md` ("What landed")
- Modify: PR #43 body

**Interfaces:**
- Consumes: all code changes from Tasks 2–5.
- Produces: a green golden gate in exact mode and M1 before/after numbers for the PR.

- [ ] **Step 1: Re-run the public-data benchmark, both estimators**

Build the current branch as the `branch` arm and the pre-Task-2 commit as the `base` arm, then run the F and H sweeps for `pydegensac` only (the field is unchanged and does not need re-measuring):

```bash
cd benchmarks
PATH="$PWD/.ab/env/bin:$PATH" CMAKE_PREFIX_PATH="$PWD/.ab/env" \
    .ab/env/bin/python -m pip install -q --no-deps --target .ab/pkg-perf ..
PYTHONPATH=$PWD/.ab/pkg-perf .ab/env/bin/python run.py f --methods pydegensac \
    --label "branch:m1-perf" --out f_m1_perf.jsonl
PYTHONPATH=$PWD/.ab/pkg-perf .ab/env/bin/python run.py h --methods pydegensac \
    --label "branch:m1-perf" --out h_m1_perf.jsonl
.ab/env/bin/python report.py f results/f_m1_perf.jsonl results/f_m1_branch.jsonl
.ab/env/bin/python report.py h results/h_m1_perf.jsonl results/h_m1_branch.jsonl
```

Expected: speed-up visible in the ms columns; no mAA difference whose paired CI excludes zero.

- [ ] **Step 2: Regenerate the golden baselines**

```bash
cd /Users/oldufo/dev/pydegensac
python3 scripts/make_golden_data.py
git status --short tests/data
```

Watch for the qualifying pair set changing (a pair failing capture sanity, as happened in `edabc68`). If a pair drops, `git rm` the stale `.npz`; if one appears, `git add` it.

- [ ] **Step 3: Review the new baselines pair by pair**

Compare old vs new inlier counts and GT error, as for `edabc68`, using the comparison against the previous commit's fixtures. The new outputs should be no worse: inlier counts flat or up, GT error flat or better. Investigate anything that moves the wrong way on both metrics at once.

- [ ] **Step 4: Verify the gate and the smoke test**

```bash
python3 -m pytest tests -q                      # expect 33 passed
cd examples && python3 simple-example.py        # expect a sane F and inlier count
```

- [ ] **Step 5: Update the write-ups and the PR**

Add a subsection to the M1 section of `docs/reports/2026-08-11-ransac-benchmark.md` with the before/after per estimator, the statistical-gate p-values, and an explicit note that the numbers are M1-only pending the maintainer's Linux re-run. Add the commits to the handoff's "What landed" table. Update the PR #43 body's speed table (via `gh api -X PATCH`, since `gh pr edit` fails on this repo with a Projects-classic error).

- [ ] **Step 6: Commit and push**

```bash
git add tests/data docs/reports
git commit -m "regenerate baselines and re-measure after the cov_mat/inlidxs work"
git push
```
