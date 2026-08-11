# Homography time-mAA

- `branch:speed-up3@0ae6b99` — h_branch.jsonl, pydegensac 0.3.0, cv2 4.13.0, poselib 2.0.5, repo e8c028a, host DM-PC
  - config: `{"EVD": {"cv2-ransac": {"px_th": 1.0, "ratio_th": null, "source": "run.py tune h (val)"}, "cv2-magsac": {"px_th": 4.0, "ratio_th": null, "source": "run.py tune h (val)"}, "poselib": {"px_th": 8.0, "ratio_th": null, "source": "run.py tune h (val)"}, "poselib-prosac": {"px_th": 4.0, "ratio_th": null, "source": "run.py tune h (val), after score-orientation fix"}, "pydegensac": {"px_th": 1.0, "ratio_th": null, "source": "run.py tune h (val)"}}, "HPatchesSeq": {"cv2-ransac": {"px_th": 16.0, "ratio_th": null, "source": "run.py tune h (val)"}, "cv2-magsac": {"px_th": 64.0, "ratio_th": null, "source": "run.py tune h (val)"}, "poselib": {"px_th": 16.0, "ratio_th": null, "source": "run.py tune h (val)"}, "poselib-prosac": {"px_th": 16.0, "ratio_th": null, "source": "run.py tune h (val), after score-orientation fix"}, "pydegensac": {"px_th": 4.0, "ratio_th": null, "source": "run.py tune h (val)"}}}`
- `base:master@08464ca` — h_base.jsonl, pydegensac 0.3.0, cv2 4.13.0, poselib 2.0.5, repo e8c028a, host DM-PC
  - config: `{"EVD": {"pydegensac": {"px_th": 1.0, "ratio_th": null, "source": "run.py tune h (val)"}}, "HPatchesSeq": {"pydegensac": {"px_th": 4.0, "ratio_th": null, "source": "run.py tune h (val)"}}}`

### EVD — best per method

| method | best mAA | 95% CI | at budget | mean ms/pair | d mAA vs poselib-prosac (paired) |
|---|---|---|---|---|---|
| poselib-prosac | 0.4500 | 0.2000-0.6875 | 100 | 1.02 | leader |
| poselib | 0.4375 | 0.1875-0.6875 | 25000 | 2.98 | -0.0126 (-0.0375, +0.0000) |
| pydegensac (branch) | 0.4250 | 0.1875-0.6625 | 1600 | 3.29 | -0.0265 (-0.1375, +0.0625) |
| cv2-magsac | 0.4000 | 0.1625-0.6375 | 6400 | 0.93 | -0.0505 (-0.1000, -0.0125) * |
| pydegensac (base) | 0.4000 | 0.1750-0.6250 | 6400 | 8.22 | -0.0512 (-0.1625, +0.0375) |
| cv2-ransac | 0.3375 | 0.0994-0.6125 | 1600 | 14.93 | -0.1131 (-0.2500, +0.0000) |

`*` = paired CI excludes zero. Each method is taken at its own best budget, which flatters every method equally.

#### pydegensac: base vs branch, same budget

| budget | mAA base | mAA branch | d mAA (95% CI, paired) | ms base | ms branch | speedup |
|---|---|---|---|---|---|---|
| 10 | 0.1750 | 0.0000 | -0.1766 (-0.4000, +0.0000) | 0.38 | 0.12 | 3.24x |
| 25 | 0.1250 | 0.1000 | -0.0252 (-0.0500, +0.0000) | 0.40 | 0.34 | 1.17x |
| 100 | 0.2125 | 0.1000 | -0.1131 (-0.2250, -0.0250) * | 0.79 | 0.50 | 1.57x |
| 400 | 0.2750 | 0.3000 | +0.0252 (+0.0000, +0.0500) | 1.31 | 1.42 | 0.92x |
| 1600 | 0.3875 | 0.4250 | +0.0369 (+0.0000, +0.0875) | 3.54 | 3.29 | 1.08x |
| 6400 | 0.4000 | 0.4000 | +0.0004 (-0.0375, +0.0375) | 8.22 | 5.74 | 1.43x |
| 25000 | 0.4000 | 0.3750 | -0.0249 (-0.0625, +0.0000) | 23.59 | 14.11 | 1.67x |
| **all** | | | | 38.23 | 25.52 | **1.50x** |

### HPatchesSeq — best per method

| method | best mAA | 95% CI | at budget | mean ms/pair | d mAA vs poselib-prosac (paired) |
|---|---|---|---|---|---|
| poselib-prosac | 0.9297 | 0.8903-0.9628 | 1600 | 7.74 | leader |
| cv2-magsac | 0.9269 | 0.8896-0.9586 | 25000 | 0.91 | -0.0027 (-0.0124, +0.0076) |
| poselib | 0.9262 | 0.8876-0.9593 | 25000 | 5.42 | -0.0034 (-0.0214, +0.0097) |
| pydegensac (branch) | 0.9172 | 0.8786-0.9503 | 25000 | 4.53 | -0.0124 (-0.0186, -0.0062) * |
| cv2-ransac | 0.9103 | 0.8710-0.9448 | 25000 | 22.65 | -0.0192 (-0.0297, -0.0090) * |
| pydegensac (base) | 0.9103 | 0.8724-0.9428 | 25000 | 6.29 | -0.0194 (-0.0290, -0.0110) * |

`*` = paired CI excludes zero. Each method is taken at its own best budget, which flatters every method equally.

#### pydegensac: base vs branch, same budget

| budget | mAA base | mAA branch | d mAA (95% CI, paired) | ms base | ms branch | speedup |
|---|---|---|---|---|---|---|
| 10 | 0.7200 | 0.7000 | -0.0196 (-0.0759, +0.0345) | 1.64 | 1.56 | 1.05x |
| 25 | 0.8041 | 0.7828 | -0.0209 (-0.0703, +0.0262) | 1.82 | 1.74 | 1.05x |
| 100 | 0.8614 | 0.8710 | +0.0091 (-0.0186, +0.0386) | 2.24 | 1.97 | 1.13x |
| 400 | 0.8883 | 0.8869 | -0.0017 (-0.0262, +0.0207) | 2.59 | 2.33 | 1.11x |
| 1600 | 0.8972 | 0.9055 | +0.0084 (-0.0014, +0.0186) | 3.03 | 2.65 | 1.14x |
| 6400 | 0.9028 | 0.9103 | +0.0077 (-0.0021, +0.0179) | 4.14 | 3.34 | 1.24x |
| 25000 | 0.9103 | 0.9172 | +0.0071 (-0.0028, +0.0172) | 6.29 | 4.53 | 1.39x |
| **all** | | | | 21.73 | 18.12 | **1.20x** |
wrote results/time_maa_h.png
