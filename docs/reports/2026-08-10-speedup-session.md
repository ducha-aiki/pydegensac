# Speed-up session report (2026-08-10, branch `speed-up3`)

Follow-up to `2026-08-10-cleanup-handoff.md`. Four commits; net result on the
golden pairs (sum of per-pair median runtimes, M-series mac, seed 42):

| estimator | before | after | speedup |
|---|---|---|---|
| `findHomography` (5 EVD pairs) | 38.4 ms | 9.2 ms | **4.2x** |
| `findFundamentalMatrix` (5 reichstag pairs) | 45.1 ms | 20.4 ms | **2.2x** |

Outputs validated at every step: bit-exact golden gate for commits 1-3,
1000-seed statistical A/B for commit 4 (which changes the sample stream and
came with a deliberate golden baseline reset).

## What was done

### 1. `resids=NULL` in bindings (`dfe0d4c`, bit-exact)

As planned in the handoff. Correct, but wall-clock-negligible at golden-pair
sizes — the residual memcpy traffic was small next to model fitting.

### 2. Profiling found the real hotspot: per-iteration RNG reseeding

`sample`-based native profiles on the golden pairs showed libc
`srandom()`+`random()` at **77% of H runtime and 51% of F runtime** — far
ahead of every candidate on the backlog (malloc churn, buffer layout were
noise). Cause: the RANSAC loops re-keyed the RNG every iteration
(`reseed_rng(rand_seed); ... rand_seed = rand()`), and each macOS `srandom()`
re-derives the 31-word BSD TYPE_3 state and discards 310 warm-up draws behind
a lock (~3.5 us per iteration).

### 3. Bit-exact local RNG (`13f68cb`, `d82348d`)

`src/pydegensac/degensac/bsd_random.c`: lock-free reimplementation of the
BSD TYPE_3 `srandom()`/`random()` pair, plus two algebraic shortcuts:

- the 310-step warm-up is linear mod 2^32, so it collapses to a precomputed
  31x31 coefficient-matrix product (961 independent multiply-adds the
  compiler vectorizes, replacing a serial dependency chain);
- the seeding LCG is Park-Miller, so `init[i] = 16807^i * seed mod 2^31-1`
  with precomputed powers and Mersenne-prime folding (serial libc-equivalent
  fallback for seed 0 or seeds with the top bit set).

Bit-exactness vs macOS libc verified over a 3M-seed sweep plus edge seeds
(0, 2^31-1, 2^31, 3e9, 2^32-1); the stream also equals glibc for every
nonzero seed. Windows keeps its existing `random`->`rand` mapping, stream
unchanged. `rand()`/`srand()` (seed chain, DegUtils) stayed libc. Reseed
cost 579 -> 62 ns; golden gate stayed green in exact mode throughout.

### 4. Remove per-iteration reseeding (`4430bc7`, behavior change + baseline reset)

The per-iteration re-key adds no entropy — results are fully determined by
the initial seed either way — so the generators are now seeded once at
estimator entry (both `exp_ranH.c` and `exp_ranF.c`; determinism via `seed`
is unchanged).

**Validation protocol**: 1000 seeded runs (seeds 1-1000) per golden pair,
before vs after, collecting per run: inlier count, median error vs GT
(reprojection px for H, symmetric epipolar px for F, over GT-consistent
correspondences at 3 px), and GT-precision of the reported inliers.
All 30 pair x metric distributions indistinguishable: two-sample KS tests,
min p = 0.108 (alpha = 0.001 with 30 tests); every mean matched within a
fraction of its std; 0 of 10,000 runs on either side failed to return a
model. Harness: scratchpad `stat_runs.py` / `compare_stats.py` (session
scratch; protocol documented here).

| pair | inliers (before / after) | GT err px (before / after) | precision (before / after) |
|---|---|---|---|
| h_adam | 110.2±12.4 / 110.8±12.4 | 1.79±0.36 / 1.78±0.35 | 0.90 / 0.90 |
| h_cat | 49.8±6.9 / 49.3±7.2 | 6.62±2.41 / 6.71±2.42 | 0.17 / 0.17 |
| h_dum | 25.8±8.2 / 26.1±7.9 | 7.54±6.23 / 7.40±5.68 | 0.55 / 0.56 |
| h_face | 82.2±13.0 / 82.6±13.9 | 1.67±0.32 / 1.67±0.32 | 0.93 / 0.93 |
| h_fox | 48.4±5.6 / 48.6±5.6 | 3.84±1.28 / 3.82±1.32 | 0.48 / 0.48 |
| f_05461164-05466646 | 63.1±7.3 / 62.8±7.4 | 0.44±0.15 / 0.45±0.15 | 0.97 / 0.96 |
| f_05534141-05545431 | 340.5±11.5 / 341.0±11.2 | 0.35±0.05 / 0.34±0.05 | 0.99 / 0.99 |
| f_05545431-05791347 | 173.0±5.2 / 172.6±5.4 | 0.45±0.07 / 0.45±0.07 | 0.99 / 0.99 |
| f_05791347-05866831 | 203.9±9.7 / 203.3±9.7 | 0.55±0.12 / 0.55±0.12 | 0.72 / 0.72 |
| f_06229406-06373813 | 189.9±7.7 / 189.9±8.1 | 0.43±0.06 / 0.43±0.07 | 0.99 / 0.99 |

**Golden baseline reset**: regenerated via `scripts/make_golden_data.py`
(self-checked determinism + GT sanity at capture). One selection change:
reichstag pair 05791347-05866831 — the borderline pair of the set (GT
precision 0.72) — no longer passes capture sanity at seeds 42/2024 under the
new stream and was replaced by 05466646-05534141, keeping 5 F pairs. Full
suite (33 tests) green in exact mode against the new baselines.

## Post-change profile (macOS, top-of-stack)

- **H**: `HDs` (error fn) is the top consumer, then `multirsampleT`,
  `nullspace`, `inlidxs`, `all_Hori_valid`. RNG cost is now minor.
- **F**: `inlidxs` (~23%), `nullspace` (~21%), `FDs` (~19%), then
  `all_ori_valid`, `rroots3`, small memmove/malloc traffic in LO.

## Remaining non-algorithmic candidates (same-results-by-construction)

Ranked by expected value from the profiles:

1. **Fuse `inlidxs` into the error loops** (F's biggest remaining item):
   error computation (`FDs`/`HDs`) writes `d[]`, then `inlidxs` re-scans it
   to collect inlier indices. Doing the threshold test in the same pass
   saves a full memory sweep per model. Same comparisons, same results.
2. **Vectorization-friendly restructuring of `FDs`/`HDs`**: per-point error
   computations are independent; keeping each point's FP operation order
   intact while helping the compiler vectorize across points is bit-exact.
   Check first whether the compiler already auto-vectorizes these loops.
3. **Hoist per-LO-iteration allocations**: `exp_inFranicustom`/`exp_inHrani`
   malloc/free inside the LO loop (visible as malloc traffic in the F
   profile); buffers can be allocated once per estimator call.
4. **Swap remaining libc `rand()`** (seed chain gone; DegUtils still draws
   `rand()` in the F degeneracy path) for a local equivalent — small, and
   only bit-exact on macOS if the macOS `rand()` LCG is replicated; likely
   not worth it (rand was ~0.2% of samples).
5. **Build flags audit**: confirm wheels compile the C core at `-O3`; no
   arch-specific flags possible for portable wheels.

Not speed, but safe and pending: bindings.cpp ~170-line dedup; delete-only
option for the dead `a<0` loop in `exp_ranF.c` (fixing the bound instead is
a behavior change needing its own baseline decision); CI pre-tag items
(cibuildwheel bump on macOS/Windows for cp314, macOS publish step).

`nullspace` (21% of F) is the main thing NOT safely optimizable: it is
hand-rolled Gaussian elimination whose FP operation order defines the golden
outputs; any reordering (or LAPACK substitution) is a numeric behavior
change. Same for `rroots3`/`slcm`.

## Cross-platform note

Local golden exactness is macOS-specific as before (CI stays in sanity
mode), but with `bsd_random.c` the *sampling* stream is now identical across
macOS and glibc Linux for nonzero seeds; residual cross-platform output
differences come from libm/FP codegen, not the RNG.
