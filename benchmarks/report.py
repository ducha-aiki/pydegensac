"""Turn the sweep jsonl into markdown tables and the time-mAA figure.

    python report.py f results/f_branch.jsonl results/f_base.jsonl
    python report.py h results/h_*.jsonl --plot results/time_maa_h.png

Every (method, budget) is one point: x = mean estimator time per pair,
y = mAA over the subset. Runs from different files are kept apart by their
``label``, so the two pydegensac builds appear as separate curves.
"""
import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

import metrics

#: Colour per curve. pydegensac's two builds are the point of the figure, so
#: they get the two ends of a warm ramp; the field is cool/neutral.
COLORS = {
    "pydegensac (base)": "#f0a35e",
    "pydegensac (branch)": "#d1462f",
    "cv2-ransac": "#1baf7a",
    "cv2-magsac": "#eda100",
    "poselib": "#2a78d6",
    "poselib-prosac": "#7a4fd1",
    # release comparison (run_releases.sh): published wheels cool, local
    # builds of the same code warm, so a wheel-vs-source gap is visible.
    "pydegensac (pypi-0.1.2)": "#9fb6c9",
    "pydegensac (pypi-0.2.1)": "#5b8db8",
    "pydegensac (pypi-0.2.2)": "#2a78d6",
    "pydegensac (local-0.2.2)": "#eda100",
    "pydegensac (local-master)": "#f0a35e",
    "pydegensac (local-branch)": "#d1462f",
}
SURFACE, PAGE = "#fcfcfb", "#f9f9f7"
INK, INK2, MUTED, GRID, BASELINE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"


def load(paths):
    """-> (rows, metas). Each row gets a ``curve`` name = method + build."""
    rows, metas = [], []
    for p in paths:
        label = None
        for line in Path(p).read_text().splitlines():
            if not line.strip():
                continue
            r = json.loads(line)
            if r.get("record") == "meta":
                label, r["file"] = r["label"], str(p)
                metas.append(r)
                continue
            r["label"] = label
            # Only pydegensac is run in more than one build; leaving the other
            # methods unsuffixed keeps them as single curves when arms merge.
            arm = (label or "").split(":")[0]
            r["curve"] = (f"pydegensac ({arm})" if r["method"] == "pydegensac"
                          else r["method"])
            rows.append(r)
    return rows, metas


#: Bootstrap resamples for the mAA confidence interval.
N_BOOT = 2000
BOOT_SEED = 0


def maa_ci(errs, problem, n_boot=N_BOOT):
    """Percentile bootstrap 95% CI of mAA over pairs.

    Backends are unseeded and the pair sets are small, so differences of a
    couple of points are not readable without this: two runs of the identical
    configuration can differ by more than the gap between two methods.
    """
    errs = np.asarray(errs, float)
    ths = metrics.F_THRESHOLDS if problem == "f" else metrics.H_THRESHOLDS
    strict = problem == "f"
    if errs.size == 0:
        return 0.0, 0.0
    idx = np.random.default_rng(BOOT_SEED).integers(
        0, errs.size, size=(n_boot, errs.size))
    sample = errs[idx]
    hits = [(sample < t).mean(1) if strict else (sample <= t).mean(1)
            for t in ths]
    boot = np.mean(hits, axis=0)
    return float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))


def _maa_over(sample, problem):
    """mAA per bootstrap row of an (n_boot, n_pairs) error sample."""
    ths = metrics.F_THRESHOLDS if problem == "f" else metrics.H_THRESHOLDS
    strict = problem == "f"
    return np.mean([(sample < t).mean(1) if strict else (sample <= t).mean(1)
                    for t in ths], axis=0)


def per_pair_errors(rows):
    """-> {(subset, curve, budget): {pair_name: error}}"""
    out = defaultdict(dict)
    for r in rows:
        err = metrics.FAIL_ERR if r["err"] is None else r["err"]
        out[(r["subset"], r["curve"], r["budget"])][r["pair"]] = err
    return out


def paired_delta_ci(errs_a, errs_b, problem, n_boot=N_BOOT):
    """95% CI of mAA(a) - mAA(b), resampling *pairs* jointly.

    Every configuration is scored on the same pairs, so the comparison is
    paired: resampling the pair set jointly cancels the shared difficulty
    that makes the marginal CIs so wide, and is what decides whether a
    difference between two configurations is real.
    """
    common = sorted(set(errs_a) & set(errs_b))
    if not common:
        return 0.0, 0.0, 0.0
    a = np.array([errs_a[k] for k in common], float)
    b = np.array([errs_b[k] for k in common], float)
    idx = np.random.default_rng(BOOT_SEED).integers(
        0, len(common), size=(n_boot, len(common)))
    delta = _maa_over(a[idx], problem) - _maa_over(b[idx], problem)
    return (float(np.percentile(delta, 2.5)),
            float(np.percentile(delta, 97.5)),
            float(np.mean(delta)))


def curves(rows, problem):
    """-> {subset: {curve: [(mean_s, median_s, maa, budget, n_fail, n,
    ci_lo, ci_hi), ...]}}"""
    maa_fn = metrics.maa_f if problem == "f" else metrics.maa_h
    groups = defaultdict(list)
    for r in rows:
        groups[(r["subset"], r["curve"], r["budget"])].append(r)

    out = defaultdict(lambda: defaultdict(list))
    for (subset, curve, budget), recs in groups.items():
        errs = [metrics.FAIL_ERR if r["err"] is None else r["err"]
                for r in recs]
        times = np.array([r["time"] for r in recs])
        lo, hi = maa_ci(errs, problem)
        out[subset][curve].append((
            float(times.mean()), float(np.median(times)), maa_fn(errs),
            budget, sum(e == metrics.FAIL_ERR for e in errs), len(recs),
            lo, hi))
    for subset in out:
        for curve in out[subset]:
            out[subset][curve].sort(key=lambda t: t[3])
    return out


def table(subset, data, unit_ms=True):
    scale = 1000.0 if unit_ms else 1.0
    lines = [f"\n### {subset}\n",
             "| method | budget | mAA | 95% CI | mean ms/pair | "
             "median ms/pair | failures |",
             "|---|---|---|---|---|---|---|"]
    for curve in sorted(data):
        for mean_t, med_t, maa, budget, n_fail, n, lo, hi in data[curve]:
            lines.append(f"| {curve} | {budget} | {maa:.4f} | "
                         f"{lo:.4f}-{hi:.4f} | {mean_t * scale:.2f} | "
                         f"{med_t * scale:.2f} | {n_fail}/{n} |")
    return "\n".join(lines)


def best_table(subset, data, per_pair, problem):
    """Peak mAA per method, each compared against the leader with a paired
    bootstrap. The marginal CI says how well the peak itself is pinned down;
    the paired column says whether the gap to the leader is real, which the
    (much wider) marginal CIs cannot answer."""
    rank = sorted(((max(pts, key=lambda t: t[2]), curve)
                   for curve, pts in data.items()),
                  key=lambda t: -t[0][2])
    leader_best, leader = rank[0]
    leader_errs = per_pair[(subset, leader, leader_best[3])]

    lines = [f"\n### {subset} — best per method\n",
             f"| method | best mAA | 95% CI | at budget | mean ms/pair | "
             f"d mAA vs {leader} (paired) |",
             "|---|---|---|---|---|---|"]
    for best, curve in rank:
        if curve == leader:
            gap = "leader"
        else:
            lo, hi, mid = paired_delta_ci(
                per_pair[(subset, curve, best[3])], leader_errs, problem)
            gap = (f"{mid:+.4f} ({lo:+.4f}, {hi:+.4f})"
                   + ("" if lo <= 0 <= hi else " *"))
        lines.append(f"| {curve} | {best[2]:.4f} | "
                     f"{best[6]:.4f}-{best[7]:.4f} | {best[3]} | "
                     f"{best[0] * 1000:.2f} | {gap} |")
    lines.append("\n`*` = paired CI excludes zero. Each method is taken at its "
                 "own best budget, which flatters every method equally.")
    return "\n".join(lines)


def speedup_table(subset, data, per_pair, problem):
    """Branch vs base pydegensac at equal iteration budget."""
    base = {p[3]: p for p in data.get("pydegensac (base)", [])}
    branch = {p[3]: p for p in data.get("pydegensac (branch)", [])}
    shared = sorted(set(base) & set(branch))
    if not shared:
        return ""
    lines = ["\n#### pydegensac: base vs branch, same budget\n",
             "| budget | mAA base | mAA branch | d mAA (95% CI, paired) | "
             "ms base | ms branch | speedup |",
             "|---|---|---|---|---|---|---|"]
    for b in shared:
        a, c = base[b], branch[b]
        lo, hi, mid = paired_delta_ci(
            per_pair[(subset, "pydegensac (branch)", b)],
            per_pair[(subset, "pydegensac (base)", b)], problem)
        sig = "" if lo <= 0 <= hi else " *"
        lines.append(f"| {b} | {a[2]:.4f} | {c[2]:.4f} | "
                     f"{mid:+.4f} ({lo:+.4f}, {hi:+.4f}){sig} | "
                     f"{a[0] * 1000:.2f} | {c[0] * 1000:.2f} | "
                     f"{a[0] / c[0]:.2f}x |")
    tot_a = sum(base[b][0] for b in shared)
    tot_c = sum(branch[b][0] for b in shared)
    lines.append(f"| **all** | | | | {tot_a * 1000:.2f} | {tot_c * 1000:.2f} | "
                 f"**{tot_a / tot_c:.2f}x** |")
    return "\n".join(lines)


def plot(all_curves, problem, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    subsets = sorted(all_curves)
    fig, axes = plt.subplots(1, len(subsets), figsize=(6.0 * len(subsets), 4.6),
                             squeeze=False, facecolor=PAGE)
    ylabel = ("mAA (1-10 deg)" if problem == "f"
              else "mAA (1-20 px, log-spaced)")
    for ax, subset in zip(axes[0], subsets):
        ax.set_facecolor(SURFACE)
        ax.set_xscale("log")
        for curve, pts in sorted(all_curves[subset].items()):
            xy = np.array([(p[0] * 1000, p[2]) for p in pts])
            color = COLORS.get(curve, MUTED)
            # Bootstrap band: without it the curves look far better separated
            # than the pair counts support.
            ax.fill_between(xy[:, 0], [p[6] for p in pts], [p[7] for p in pts],
                            color=color, alpha=0.12, linewidth=0, zorder=2)
            ax.plot(xy[:, 0], xy[:, 1], "-o", color=color,
                    linewidth=2, markersize=6, label=curve, zorder=3,
                    markeredgecolor=SURFACE, markeredgewidth=1.2)
        ax.set_title(subset, color=INK, fontsize=10.5)
        ax.set_xlabel("mean time per pair (ms, log scale)", color=MUTED,
                      fontsize=9)
        ax.grid(True, color=GRID, linewidth=0.75, zorder=0)
        ax.tick_params(colors=MUTED, labelsize=8.5)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(BASELINE)
        ax.margins(x=0.12)
    axes[0][0].set_ylabel(ylabel, color=MUTED, fontsize=9)
    # One legend below the panels: curves converge in the lower right, which
    # is exactly where an in-axes legend would sit.
    handles, labels = axes[0][0].get_legend_handles_labels()
    # Give every legend column ~2.7in; release comparisons have six long
    # labels that overrun a fixed column count on a single-panel figure.
    ncol = max(1, min(len(labels), int(6.0 * len(subsets) / 2.7)))
    nrow = -(-len(labels) // ncol)
    fig.legend(handles, labels, loc="lower center", ncol=ncol,
               fontsize=8.5, frameon=False, labelcolor=INK2,
               bbox_to_anchor=(0.5, 0.0))
    fig.suptitle(("Fundamental matrix" if problem == "f" else "Homography")
                 + ": accuracy vs. compute budget", color=INK, fontsize=12.5,
                 y=0.985)
    # The bands are marginal CIs — they show how loosely each curve is pinned
    # down, but overlapping bands do NOT mean two methods are tied. That is a
    # paired question, answered in the tables.
    fig.text(0.5, 0.925, "bands: marginal 95% bootstrap CI; method gaps are "
             "tested paired (tables)", ha="center", fontsize=8, color=MUTED)
    fig.tight_layout(rect=(0, 0.055 * nrow, 1, 0.905))
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    print(f"wrote {out_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("problem", choices=["f", "h"])
    ap.add_argument("results", nargs="+")
    ap.add_argument("--plot", help="write the time-mAA figure here")
    ap.add_argument("--full", action="store_true",
                    help="also print the per-budget table")
    args = ap.parse_args()

    rows, metas = load(args.results)
    all_curves = curves(rows, args.problem)
    per_pair = per_pair_errors(rows)

    print(f"# {'Fundamental' if args.problem == 'f' else 'Homography'} "
          f"time-mAA\n")
    for m in metas:
        print(f"- `{m['label']}` — {Path(m['file']).name}, "
              f"pydegensac {m['versions'].get('pydegensac')}, "
              f"cv2 {m['versions'].get('cv2')}, "
              f"poselib {m['versions'].get('poselib')}, "
              f"repo {m['git_sha']}, host {m['host']}")
        if "config" in m:
            print(f"  - config: `{json.dumps(m['config'])}`")

    for subset in sorted(all_curves):
        print(best_table(subset, all_curves[subset], per_pair, args.problem))
        print(speedup_table(subset, all_curves[subset], per_pair, args.problem))
        if args.full:
            print(table(subset, all_curves[subset]))

    if args.plot:
        plot(all_curves, args.problem, args.plot)


if __name__ == "__main__":
    main()
