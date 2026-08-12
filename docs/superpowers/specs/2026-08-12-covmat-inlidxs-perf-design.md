# `cov_mat` and `inlidxs` — design (2026-08-12)

## Purpose

The `speed-up3` profile that guided the RNG work was taken on a build where
every LAPACK call was preprocessed away (`f349a6c`), so local optimisation was
not running. With LO restored the profile is different, and the two functions
it now points at are both structurally slow rather than merely hot.

Self time, M1, fixed build, `sample` over ~15 s of estimation on real pairs:

| F (10k iters, st_peters pairs) | | H (25k iters, HPatchesSeq) | |
|---|---|---|---|
| `inlidxs` | ~46% | `cov_mat` | ~60% |
| `FDs` | ~18% | `HDs` | ~15% |
| `cov_mat` | ~8% | `inlidxs` | ~6% |
| LAPACK (DSTEQR/DLASR/…) | ~8% | `normu` | ~5% |

The error functions (`FDs`/`HDs`) are **not** the target: the 2026-08-10
addendum measured them as already auto-vectorised at width 2, division-limited,
and this profile agrees they are second-order.

Both changes preserve the algorithm exactly — same models, same iteration
sequence, same decisions — and change only floating-point rounding.

## Non-goals

- No algorithmic change: no adaptive LO, no sampling change, no threshold
  change.
- No SIMD intrinsics and no architecture-specific code. Wheels are portable.
- Not bit-exactness. This is a deliberate numerics change, gated statistically
  (see below), following the precedent set by the per-iteration reseed removal
  (`4430bc7`).

## Change 1 — `cov_mat`: one pass, not 45

`utools.c:170`. Computes `Cv = Zᵀ Z` for `Z` of shape `len x siz`, row-major:

```c
for (i=0; i<siz; i++)
  for (j=0; j<=i; j++) {
      val = 0;
      for (k=0; k<lenM; k+=siz) val += Z[k+i] * Z[k+j];   /* serial chain */
      Cv[siz*i + j] = val; Cv[i + siz*j] = val;
  }
```

Two problems, both structural:

1. **45 passes over the data** (for `siz=9`), each touching 2 of every 9
   doubles.
2. **Each pass is one serial FP dependency chain.** At ~4-cycle FMA latency the
   inner loop is latency-bound, not throughput-bound: ~`4*len` cycles per pass,
   ~`180*len` cycles total.

The hot call site is `exp_ranH.c`'s `cov_mat(C, A, 10, 9)` — inside the
per-iteration MCE solve, so it runs on **every** RANSAC iteration. At `len=10`
that is ~1800 cycles per call against ~112 if the 45 accumulators are
independent.

**Approach: single pass over points, accumulating all `siz*(siz+1)/2` unique
entries.** The accumulators are mutually independent, so the FMA units
saturate; `Z` is read once.

`dsyrk` was considered and rejected as the primary form: the hot sites have
`len` of 8–10, where BLAS call overhead is comparable to the entire
computation. It remains a candidate for the large-`len` sites (`u2f`/`u2h`,
over all inliers) and will be added behind a length threshold **only if
measured faster there**.

Constraints:

- `siz` is not always 9 — `exp_ranH.c:751` passes `nullsize`. The
  implementation must stay correct for arbitrary `siz`, with no fixed-size
  assumption.
- Both triangles of `Cv` are filled, as now.

Numeric effect: summation order over points is unchanged (still ascending `k`);
what changes is that intermediates live in different registers and the compiler
is free to contract differently. Rounding may differ in the last bits.

## Change 2 — `inlidxs` / `truncQuad`

`rtools.c:163` and `rtools.c:231`. Per point, the loop makes a cross-TU call to
`truncQuad`, which recomputes `thr*9/4`, branches twice, and **divides**:

```c
for (i = 0; i < len; ++i) {
    s.J += truncQuad(err[i], th);          /* call + 2 branches + divide */
    if (err[i] <= th) { inl[s.I] = i; ++(s.I); }
}
```

Three steps, deliberately separable so the gate can attribute any movement:

1. `truncQuad` becomes `static inline` in `rtools.h`. **Bit-exact.**
2. Branchless compress store: `inl[s.I] = i; s.I += (err[i] <= th);`.
   **Bit-exact** — no FP operation changes.
3. Hoisted reciprocal: `inv = 1.0/(th*9/4)` once per call, then
   `J += fmax(0.0, 1.0 - err[i]*inv)`. **This is the numerics change.**

Step 3 is equivalent in exact arithmetic: `truncQuad` returns 0 exactly when
`epsilon >= thr*9/4`, which is where `1 - epsilon*inv <= 0`. The `thr == 0`
guard moves out of the loop.

`inl` must have room for `len` entries for the branchless store to be safe —
true at every call site today (buffers are allocated at `len`), and to be
asserted rather than assumed.

## Validation — `scripts/stat_ab.py` (new, committed)

The 2026-08-10 protocol was ad hoc and was not kept. It is rebuilt here as a
committed script, because this is now the second deliberate numerics change and
there will be others.

Per golden pair, 1000 seeded runs (seeds 1–1000) against a reference build,
collecting per run:

- inlier count,
- median error vs ground truth (reprojection px for H, symmetric epipolar px
  for F, over GT-consistent correspondences at 3 px),
- GT precision of the reported inliers.

Two-sample KS test per (pair × metric) — 30 comparisons over the 10 golden
pairs — reporting the minimum p-value. The 2026-08-10 run's min p was 0.108;
anything comparable passes, a small p is a stop-and-investigate.

The script takes two importable builds (via `PYTHONPATH`, as `benchmarks/`
already does) so it can compare any two commits, not just this change.

## Acceptance

1. `scripts/stat_ab.py` shows no distributional difference on either estimator.
2. Public-data benchmark re-run (F and H, paired bootstrap): no mAA difference
   whose CI excludes zero.
3. Measured speed-up reported per estimator, M1.
4. Golden baselines regenerated, reviewed pair by pair (inlier counts and GT
   error compared old vs new, as for `edabc68`).
5. Example smoke test passes.

## Risks

- **Only macOS/M1 is measurable in this session.** The structural problem is
  platform-independent and x86-64 should benefit at least as much, but the
  reported numbers will be M1 until the maintainer re-runs on Linux.
- Second baseline regeneration inside PR #43. Both are deliberate and
  documented; a reviewer sees two resets in one branch.
- If the single-pass `cov_mat` spills registers at large `siz`, it could
  regress the large-`len` sites. Measured per call site, not assumed.
