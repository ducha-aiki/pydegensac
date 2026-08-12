"""Freeze the benchmark correspondences into a numpy-only .npz.

``iso_bench.py`` must not import cv2/h5py: they dominate the process and their
own builds differ between environments, which is exactly what is being
measured. Run once, from anywhere.
"""
import sys
from pathlib import Path

import numpy as np

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH))
import data  # noqa: E402

OUT = BENCH / ".ab" / "iso" / "pairs.npz"
out = {}

# H: HPatchesSeq test. pydegensac's tuned px_th is 4.0 with no ratio filter,
# so every pair goes in whole.
n = 0
for p in data.iter_pairs_h("HPatchesSeq", "test"):
    out[f"h/{n}/1"] = p["pts1"]
    out[f"h/{n}/2"] = p["pts2"]
    n += 1
print(f"h: {n} pairs")

# F: st_peters_square. The tuned ratio_th of 0.8 is applied here rather than in
# the driver, so the driver stays pure timing.
m = 0
for p in data.iter_pairs_f():
    keep = p["scores"] <= 0.8
    if keep.sum() < 8:
        continue
    out[f"f/{m}/1"] = np.ascontiguousarray(p["pts1"][keep])
    out[f"f/{m}/2"] = np.ascontiguousarray(p["pts2"][keep])
    m += 1
    if m == 100:  # 100 pairs is plenty for a timing loop
        break
print(f"f: {m} pairs")

OUT.parent.mkdir(parents=True, exist_ok=True)
np.savez(OUT, **out)
print(f"wrote {OUT}")
