"""Median ms/pair per arm, relative to whichever arm was listed first.

    python summarize.py h.jsonl f.jsonl
"""
import collections
import json
import statistics
import sys

for path in sys.argv[1:]:
    rows = [json.loads(l) for l in open(path) if l.startswith("{")]
    by = collections.defaultdict(list)
    checks = collections.defaultdict(set)
    for r in rows:
        by[r["label"]].append(r["ms_per_pair"])
        checks[r["label"]].add(r["inliers"])
    base = statistics.median(next(iter(by.values())))
    print(f"== {path}")
    print(f"{'arm':22s} {'n':>2} {'median':>9} {'min':>8} {'max':>8}  vs-first")
    for k, v in by.items():
        m = statistics.median(v)
        print(f"{k:22s} {len(v):2d} {m:9.3f} {min(v):8.3f} {max(v):8.3f}  "
              f"x{m / base:.3f}")
    # If the arms disagree here they are not doing the same work, and none of
    # the timings above mean anything.
    allsums = set().union(*checks.values())
    print(f"inlier checksum: {'IDENTICAL ' if len(allsums) == 1 else 'DIFFERS '}"
          f"{sorted(allsums)}\n")
