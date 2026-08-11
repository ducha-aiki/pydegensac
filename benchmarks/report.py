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


def curves(rows, problem):
    """-> {subset: {curve: [(mean_s, median_s, maa, budget, n_fail), ...]}}"""
    maa_fn = metrics.maa_f if problem == "f" else metrics.maa_h
    groups = defaultdict(list)
    for r in rows:
        groups[(r["subset"], r["curve"], r["budget"])].append(r)

    out = defaultdict(lambda: defaultdict(list))
    for (subset, curve, budget), recs in groups.items():
        errs = [metrics.FAIL_ERR if r["err"] is None else r["err"]
                for r in recs]
        times = np.array([r["time"] for r in recs])
        out[subset][curve].append((
            float(times.mean()), float(np.median(times)), maa_fn(errs),
            budget, sum(e == metrics.FAIL_ERR for e in errs), len(recs)))
    for subset in out:
        for curve in out[subset]:
            out[subset][curve].sort(key=lambda t: t[3])
    return out


def table(subset, data, unit_ms=True):
    scale = 1000.0 if unit_ms else 1.0
    lines = [f"\n### {subset}\n",
             "| method | budget | mAA | mean ms/pair | median ms/pair | failures |",
             "|---|---|---|---|---|---|"]
    for curve in sorted(data):
        for mean_t, med_t, maa, budget, n_fail, n in data[curve]:
            lines.append(f"| {curve} | {budget} | {maa:.4f} | "
                         f"{mean_t * scale:.2f} | {med_t * scale:.2f} | "
                         f"{n_fail}/{n} |")
    return "\n".join(lines)


def best_table(subset, data):
    """Peak mAA per curve, and the cheapest budget within 0.002 mAA of it."""
    lines = [f"\n### {subset} — best per method\n",
             "| method | best mAA | at budget | mean ms/pair | "
             "cheapest within 0.002 mAA | ms there |",
             "|---|---|---|---|---|---|"]
    rank = []
    for curve, pts in data.items():
        best = max(pts, key=lambda t: t[2])
        near = min((p for p in pts if p[2] >= best[2] - 0.002),
                   key=lambda t: t[0])
        rank.append((best[2], curve, best, near))
    for maa, curve, best, near in sorted(rank, reverse=True):
        lines.append(f"| {curve} | {maa:.4f} | {best[3]} | "
                     f"{best[0] * 1000:.2f} | {near[3]} | "
                     f"{near[0] * 1000:.2f} |")
    return "\n".join(lines)


def speedup_table(data):
    """Branch vs base pydegensac at equal iteration budget."""
    base = {p[3]: p for p in data.get("pydegensac (base)", [])}
    branch = {p[3]: p for p in data.get("pydegensac (branch)", [])}
    shared = sorted(set(base) & set(branch))
    if not shared:
        return ""
    lines = ["\n#### pydegensac: base vs branch, same budget\n",
             "| budget | mAA base | mAA branch | ms base | ms branch | speedup |",
             "|---|---|---|---|---|---|"]
    for b in shared:
        a, c = base[b], branch[b]
        lines.append(f"| {b} | {a[2]:.4f} | {c[2]:.4f} | {a[0] * 1000:.2f} | "
                     f"{c[0] * 1000:.2f} | {a[0] / c[0]:.2f}x |")
    tot_a = sum(base[b][0] for b in shared)
    tot_c = sum(branch[b][0] for b in shared)
    lines.append(f"| **all** | | | {tot_a * 1000:.2f} | {tot_c * 1000:.2f} | "
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
            ax.plot(xy[:, 0], xy[:, 1], "-o", color=COLORS.get(curve, MUTED),
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
    fig.legend(handles, labels, loc="lower center", ncol=len(labels),
               fontsize=8.5, frameon=False, labelcolor=INK2,
               bbox_to_anchor=(0.5, 0.0))
    fig.suptitle(("Fundamental matrix" if problem == "f" else "Homography")
                 + ": accuracy vs. compute budget", color=INK, fontsize=12.5)
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
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
        print(best_table(subset, all_curves[subset]))
        print(speedup_table(all_curves[subset]))
        if args.full:
            print(table(subset, all_curves[subset]))

    if args.plot:
        plot(all_curves, args.problem, args.plot)


if __name__ == "__main__":
    main()
