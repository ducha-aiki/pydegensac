# Handoff: post-cleanup state and speed-up roadmap (2026-08-10)

Audience: the next working session on pydegensac. Everything below is on `master` after PR #42 (merge commit `7760039`); version bumped to 0.3.0.

## What just landed (PR #42)

- **Golden regression safety net.** `tests/test_golden_regression.py` + committed baselines in `tests/data/`: 5 EVD homography pairs (inputs are MODS pre-extracted tentatives from `EVD_tentatives.zip`, NOT SIFT — vanilla SIFT fails on EVD) and 5 IMC reichstag F pairs (SIFT matches via imc-2021-simple's extractor), each with GT (H, or F + K/R/T), captured at seeds 42 and 2024.
- **Dead code removed, bit-equivalent** (~5,900 lines): `ranF.c/h`, `ranH2el.c/h`, most of `ranH.c` (only `inHrani`+`iterH` survive — used by `DegUtils.c`'s degeneracy check), dead duplicates in `exp_ranH.c`/`exp_ranF.c` (`exp_iterH`, `exp_inHrani`, `hMCEs`, `exp_ransacF`, `exp_ransacFcustom`, `exp_inFrani`, `exp_iterF`), matutls trimmed 61→11 files. Live chain: `bindings.cpp` → `exp_ransacHcustomLAF` / `exp_ransacFcustomLAF` (+ `*custom` helpers, `HcloseToSingular`, hash table).
- **Audit bugs fixed** (TDD, `tests/test_input_validation.py`): bindings zero-init H/F outputs, `c_style|forcecast`, ndim/width/[N,2]-vs-[N,6] checks, LAF-with-Nx2 guard, malloc checks, `PYBIND11_MODULE`; utils.py failure-mask type + typed excepts; setup.py py3.12-safe; CMake honors `--debug`; `dgeqp3_` declared in `lapwrap.h` (GCC 14 hard-errors on implicit declarations).
- **CI runs the test suite** in all four cibuildwheel jobs (sanity mode). Full matrix green: 20 wheels × 33 tests. Linux is split: manylinux2014 (cp38-cp313) + manylinux_2_28 (cp314, cibuildwheel 3.4.1 — the old 2.21.2 pin silently never built cp314 anywhere).

## The iron rule: the golden gate

Any change to C code, bindings, or build flags must keep the golden tests bit-exact:

```bash
rm -rf build/ && pip install . --force-reinstall --no-deps -q && python -m pytest tests -q
```

- `rm -rf build/` matters: pip reuses stale object files and you can "pass" against an old binary.
- Local default is exact mode. CI sets `PYDEGENSAC_GOLDEN_EXACT=0` (determinism + inlier-floor + GT-agreement) because `rand()`/FP codegen differ per libc — do not expect cross-platform bit-exactness.
- Never loosen a golden test or regenerate baselines to make a change pass. Regeneration (`scripts/make_golden_data.py`) is a deliberate baseline reset: it only runs on this machine (hardcoded `/Users/oldufo/dev/imc-2021-simple` path, cached datasets in gitignored `.cache/`), needs the local build installed first, and self-checks determinism + GT sanity at capture.

## Speed-up backlog (the next planned work), in suggested order

1. **`resids=NULL` in bindings** — verified pure overhead: both binding call sites pass `&resids`, making the C core realloc + memcpy full-length residual arrays every LO iteration into a buffer the bindings free unread (`bindings.cpp` frees without reading). The live `*customLAF` functions already NULL-guard. Expected: allocation/memcpy traffic gone, outputs bit-identical → gate must stay green in exact mode.
2. **Profile before further work** — after (1), profile `findFundamentalMatrix`/`findHomography` on the golden pairs (plus a larger synthetic set); candidates: the per-iteration `malloc`/`free` churn, `data_out` allocation, error-buffer layout.
3. **`exp_ranF.c` ~line 796**: `for (a=0;a<0;a++) H_best[a] = Hbest[a];` — never-executing loop; the H-from-degeneracy output (`HinF` in bindings, uninitialized but dead) is never produced. Either fix the bound (BEHAVIOR CHANGE — golden masks/models could change; would need a deliberate baseline decision) or delete the dead plumbing (bit-equivalent).
4. **The `(IDEA: ...)` note at `exp_ranF.h:10`** — adaptive LO-triggering heuristic, never implemented. Owner explicitly wants it evaluated. This is an algorithmic change: benchmark accuracy (IMC pairs, EVD) + speed against current behavior; not gated on bit-equivalence, gated on measured quality.
5. **Bindings dedup** — ~170 duplicated lines between the two binding functions; safe to unify once (1) settles, gate stays exact.

## Before tagging v0.3.0 (wheel release)

- Bump `cibuildwheel` in the **macOS and Windows** jobs (still 2.21.2 → they silently skip cp314; Linux already fixed). Watch for platform drift on first run — every CI issue so far was ecosystem drift, not our code.
- Add the missing **macOS release-publish step** (linux/windows publish wheels on tags, macOS never did — pre-existing gap).
- Optional hygiene: `concurrency:` group in the workflow; the `numpy<2.3` / `opencv-python-headless==4.*` CI test pins are environment-shaped (documented inline in the workflow) — revisit when manylinux2014 retires.

## Pitfalls learned this session

- cibuildwheel installs the wheel (and its `install_requires`) BEFORE `CIBW_TEST_REQUIRES` — constraints for install-time deps must go in `CIBW_BEFORE_TEST`.
- cmd.exe eats `<` in requirement specs — use metacharacter-free pins (`==4.*`) in CIBW variables.
- `imc2021.download.download_dataset` has no request timeout (can hang mid-download) and breaks on py3.11 (`extractall(filter=)`); `scripts/make_golden_data.py` works around both, but only when the tarball is already cached.
- Old cibuildwheel silently drops unknown build identifiers — a "green" matrix can be quietly not building what you think.

## Where things are

- Plan (13 tasks, all complete): `docs/superpowers/plans/2026-08-10-big-cleanup.md`
- Audit summary + conventions: `CLAUDE.md` (architecture section reflects the post-cleanup file set)
- This session's task reports/reviews lived in `.superpowers/sdd/2026-08-10-big-cleanup/` (git-ignored scratch; deleted after merge — git history is the record)
- Downstream consumer to keep working: `/Users/oldufo/dev/imc-2021-simple` — uses only `findFundamentalMatrix(Nx2, px_th=, conf=0.9999, max_iters=)` and needs `mask.ravel().astype(bool)` to work
