# Homography time-mAA

- `branch:m1-f349a6c` — h_m1_branch.jsonl, pydegensac 0.3.0, cv2 4.13.0, poselib 2.0.5, repo f349a6c, host Dmytros-MacBook-Air.local
  - config: `{"EVD": {"cv2-ransac": {"px_th": 1.0, "ratio_th": null, "source": "run.py tune h (val)"}, "cv2-magsac": {"px_th": 4.0, "ratio_th": null, "source": "run.py tune h (val)"}, "poselib": {"px_th": 8.0, "ratio_th": null, "source": "run.py tune h (val)"}, "poselib-prosac": {"px_th": 4.0, "ratio_th": null, "source": "run.py tune h (val), after score-orientation fix"}, "pydegensac": {"px_th": 1.0, "ratio_th": null, "source": "run.py tune h (val)"}}, "HPatchesSeq": {"cv2-ransac": {"px_th": 16.0, "ratio_th": null, "source": "run.py tune h (val)"}, "cv2-magsac": {"px_th": 64.0, "ratio_th": null, "source": "run.py tune h (val)"}, "poselib": {"px_th": 16.0, "ratio_th": null, "source": "run.py tune h (val)"}, "poselib-prosac": {"px_th": 16.0, "ratio_th": null, "source": "run.py tune h (val), after score-orientation fix"}, "pydegensac": {"px_th": 4.0, "ratio_th": null, "source": "run.py tune h (val)"}}}`
- `base:m1-master+lapackfix` — h_m1_basefix.jsonl, pydegensac 0.3.0, cv2 4.13.0, poselib 2.0.5, repo f349a6c, host Dmytros-MacBook-Air.local
  - config: `{"EVD": {"pydegensac": {"px_th": 1.0, "ratio_th": null, "source": "run.py tune h (val)"}}, "HPatchesSeq": {"pydegensac": {"px_th": 4.0, "ratio_th": null, "source": "run.py tune h (val)"}}}`

### EVD — best per method

| method | best mAA | 95% CI | at budget | mean ms/pair | d mAA vs poselib-prosac (paired) |
|---|---|---|---|---|---|
| poselib-prosac | 0.4500 | 0.2000-0.6875 | 100 | 0.82 | leader |
| poselib | 0.4375 | 0.1875-0.6875 | 25000 | 2.90 | -0.0126 (-0.0375, +0.0000) |
| pydegensac (base) | 0.4250 | 0.1875-0.6750 | 6400 | 9.65 | -0.0259 (-0.1750, +0.1125) |
| cv2-magsac | 0.4125 | 0.1750-0.6500 | 6400 | 1.03 | -0.0379 (-0.0750, -0.0125) * |
| pydegensac (branch) | 0.3625 | 0.1625-0.5625 | 6400 | 4.09 | -0.0884 (-0.2000, +0.0250) |
| cv2-ransac | 0.3375 | 0.0994-0.6125 | 1600 | 17.79 | -0.1131 (-0.2500, +0.0000) |

`*` = paired CI excludes zero. Each method is taken at its own best budget, which flatters every method equally.

#### pydegensac: base vs branch, same budget

| budget | mAA base | mAA branch | d mAA (95% CI, paired) | ms base | ms branch | speedup |
|---|---|---|---|---|---|---|
| 10 | 0.1750 | 0.0625 | -0.1135 (-0.2875, -0.0125) * | 0.47 | 0.33 | 1.41x |
| 25 | 0.2250 | 0.0875 | -0.1387 (-0.3375, -0.0125) * | 0.41 | 0.25 | 1.60x |
| 100 | 0.2500 | 0.2625 | +0.0125 (-0.0750, +0.1125) | 0.83 | 0.58 | 1.42x |
| 400 | 0.3500 | 0.2875 | -0.0625 (-0.1875, +0.0500) | 1.99 | 0.87 | 2.29x |
| 1600 | 0.3625 | 0.3500 | -0.0149 (-0.1375, +0.1125) | 3.65 | 1.98 | 1.84x |
| 6400 | 0.4250 | 0.3625 | -0.0625 (-0.1750, +0.0500) | 9.65 | 4.09 | 2.36x |
| 25000 | 0.4250 | 0.3625 | -0.0625 (-0.1750, +0.0500) | 29.93 | 9.47 | 3.16x |
| **all** | | | | 46.93 | 17.59 | **2.67x** |

### HPatchesSeq — best per method

| method | best mAA | 95% CI | at budget | mean ms/pair | d mAA vs poselib-prosac (paired) |
|---|---|---|---|---|---|
| poselib-prosac | 0.9297 | 0.8903-0.9628 | 1600 | 7.80 | leader |
| cv2-magsac | 0.9269 | 0.8896-0.9586 | 25000 | 0.95 | -0.0027 (-0.0124, +0.0076) |
| poselib | 0.9262 | 0.8876-0.9593 | 25000 | 5.19 | -0.0034 (-0.0214, +0.0097) |
| pydegensac (base) | 0.9124 | 0.8745-0.9462 | 25000 | 8.23 | -0.0171 (-0.0276, -0.0076) * |
| cv2-ransac | 0.9103 | 0.8710-0.9448 | 25000 | 28.00 | -0.0192 (-0.0297, -0.0090) * |
| pydegensac (branch) | 0.9062 | 0.8648-0.9421 | 25000 | 5.23 | -0.0234 (-0.0414, -0.0103) * |

`*` = paired CI excludes zero. Each method is taken at its own best budget, which flatters every method equally.

#### pydegensac: base vs branch, same budget

| budget | mAA base | mAA branch | d mAA (95% CI, paired) | ms base | ms branch | speedup |
|---|---|---|---|---|---|---|
| 10 | 0.6966 | 0.7276 | +0.0322 (-0.0138, +0.0814) | 2.75 | 2.73 | 1.01x |
| 25 | 0.7814 | 0.7890 | +0.0085 (-0.0331, +0.0524) | 2.88 | 2.80 | 1.03x |
| 100 | 0.8386 | 0.8338 | -0.0052 (-0.0414, +0.0324) | 3.20 | 3.12 | 1.02x |
| 400 | 0.8772 | 0.8828 | +0.0056 (-0.0145, +0.0290) | 3.77 | 3.63 | 1.04x |
| 1600 | 0.8876 | 0.9000 | +0.0126 (-0.0069, +0.0352) | 4.23 | 3.91 | 1.08x |
| 6400 | 0.9117 | 0.9055 | -0.0062 (-0.0262, +0.0110) | 5.51 | 4.41 | 1.25x |
| 25000 | 0.9124 | 0.9062 | -0.0063 (-0.0262, +0.0103) | 8.23 | 5.23 | 1.57x |
| **all** | | | | 30.57 | 25.83 | **1.18x** |
wrote results/time_maa_h_m1_equalcorrectness.png
