# Homography time-mAA

- `branch:m1-f349a6c` — h_m1_branch.jsonl, pydegensac 0.3.0, cv2 4.13.0, poselib 2.0.5, repo f349a6c, host Dmytros-MacBook-Air.local
  - config: `{"EVD": {"cv2-ransac": {"px_th": 1.0, "ratio_th": null, "source": "run.py tune h (val)"}, "cv2-magsac": {"px_th": 4.0, "ratio_th": null, "source": "run.py tune h (val)"}, "poselib": {"px_th": 8.0, "ratio_th": null, "source": "run.py tune h (val)"}, "poselib-prosac": {"px_th": 4.0, "ratio_th": null, "source": "run.py tune h (val), after score-orientation fix"}, "pydegensac": {"px_th": 1.0, "ratio_th": null, "source": "run.py tune h (val)"}}, "HPatchesSeq": {"cv2-ransac": {"px_th": 16.0, "ratio_th": null, "source": "run.py tune h (val)"}, "cv2-magsac": {"px_th": 64.0, "ratio_th": null, "source": "run.py tune h (val)"}, "poselib": {"px_th": 16.0, "ratio_th": null, "source": "run.py tune h (val)"}, "poselib-prosac": {"px_th": 16.0, "ratio_th": null, "source": "run.py tune h (val), after score-orientation fix"}, "pydegensac": {"px_th": 4.0, "ratio_th": null, "source": "run.py tune h (val)"}}}`
- `base:m1-master` — h_m1_base.jsonl, pydegensac 0.3.0, cv2 4.13.0, poselib 2.0.5, repo f349a6c, host Dmytros-MacBook-Air.local
  - config: `{"EVD": {"pydegensac": {"px_th": 1.0, "ratio_th": null, "source": "run.py tune h (val)"}}, "HPatchesSeq": {"pydegensac": {"px_th": 4.0, "ratio_th": null, "source": "run.py tune h (val)"}}}`

### EVD — best per method

| method | best mAA | 95% CI | at budget | mean ms/pair | d mAA vs poselib-prosac (paired) |
|---|---|---|---|---|---|
| poselib-prosac | 0.4500 | 0.2000-0.6875 | 100 | 0.82 | leader |
| poselib | 0.4375 | 0.1875-0.6875 | 25000 | 2.90 | -0.0126 (-0.0375, +0.0000) |
| cv2-magsac | 0.4125 | 0.1750-0.6500 | 6400 | 1.03 | -0.0379 (-0.0750, -0.0125) * |
| pydegensac (branch) | 0.3625 | 0.1625-0.5625 | 6400 | 4.09 | -0.0884 (-0.2000, +0.0250) |
| cv2-ransac | 0.3375 | 0.0994-0.6125 | 1600 | 17.79 | -0.1131 (-0.2500, +0.0000) |
| pydegensac (base) | 0.3250 | 0.1250-0.5625 | 6400 | 8.12 | -0.1253 (-0.2750, +0.0250) |

`*` = paired CI excludes zero. Each method is taken at its own best budget, which flatters every method equally.

#### pydegensac: base vs branch, same budget

| budget | mAA base | mAA branch | d mAA (95% CI, paired) | ms base | ms branch | speedup |
|---|---|---|---|---|---|---|
| 10 | 0.0375 | 0.0625 | +0.0253 (+0.0000, +0.0750) | 0.17 | 0.33 | 0.52x |
| 25 | 0.0875 | 0.0875 | +0.0001 (-0.1500, +0.1500) | 0.15 | 0.25 | 0.57x |
| 100 | 0.1625 | 0.2625 | +0.1007 (+0.0000, +0.2250) | 0.30 | 0.58 | 0.51x |
| 400 | 0.1375 | 0.2875 | +0.1508 (+0.0375, +0.2750) * | 0.79 | 0.87 | 0.91x |
| 1600 | 0.1875 | 0.3500 | +0.1607 (+0.0500, +0.3000) * | 2.56 | 1.98 | 1.29x |
| 6400 | 0.3250 | 0.3625 | +0.0370 (-0.0500, +0.1125) | 8.12 | 4.09 | 1.98x |
| 25000 | 0.3125 | 0.3625 | +0.0491 (-0.0500, +0.1500) | 28.06 | 9.47 | 2.96x |
| **all** | | | | 40.14 | 17.59 | **2.28x** |

### HPatchesSeq — best per method

| method | best mAA | 95% CI | at budget | mean ms/pair | d mAA vs poselib-prosac (paired) |
|---|---|---|---|---|---|
| poselib-prosac | 0.9297 | 0.8903-0.9628 | 1600 | 7.80 | leader |
| cv2-magsac | 0.9269 | 0.8896-0.9586 | 25000 | 0.95 | -0.0027 (-0.0124, +0.0076) |
| poselib | 0.9262 | 0.8876-0.9593 | 25000 | 5.19 | -0.0034 (-0.0214, +0.0097) |
| cv2-ransac | 0.9103 | 0.8710-0.9448 | 25000 | 28.00 | -0.0192 (-0.0297, -0.0090) * |
| pydegensac (branch) | 0.9062 | 0.8648-0.9421 | 25000 | 5.23 | -0.0234 (-0.0414, -0.0103) * |
| pydegensac (base) | 0.6924 | 0.6503-0.7317 | 25000 | 5.79 | -0.2375 (-0.2690, -0.2069) * |

`*` = paired CI excludes zero. Each method is taken at its own best budget, which flatters every method equally.

#### pydegensac: base vs branch, same budget

| budget | mAA base | mAA branch | d mAA (95% CI, paired) | ms base | ms branch | speedup |
|---|---|---|---|---|---|---|
| 10 | 0.3807 | 0.7276 | +0.3483 (+0.2862, +0.4124) * | 0.56 | 2.73 | 0.21x |
| 25 | 0.5352 | 0.7890 | +0.2549 (+0.2007, +0.3124) * | 0.60 | 2.80 | 0.21x |
| 100 | 0.5924 | 0.8338 | +0.2417 (+0.2007, +0.2848) * | 0.75 | 3.12 | 0.24x |
| 400 | 0.6290 | 0.8828 | +0.2540 (+0.2131, +0.2966) * | 1.00 | 3.63 | 0.28x |
| 1600 | 0.6566 | 0.9000 | +0.2439 (+0.2034, +0.2869) * | 1.45 | 3.91 | 0.37x |
| 6400 | 0.6807 | 0.9055 | +0.2252 (+0.1917, +0.2593) * | 2.65 | 4.41 | 0.60x |
| 25000 | 0.6924 | 0.9062 | +0.2141 (+0.1834, +0.2462) * | 5.79 | 5.23 | 1.11x |
| **all** | | | | 12.80 | 25.83 | **0.50x** |
wrote results/time_maa_h_m1.png
