# Homography time-mAA

- `branch:m1-opt@15445ef` — h_m1_opt.jsonl, pydegensac 0.3.0, cv2 4.13.0, poselib 2.0.5, repo aa111c8, host Dmytros-MacBook-Air.local
  - config: `{"EVD": {"cv2-ransac": {"px_th": 1.0, "ratio_th": null, "source": "run.py tune h (val)"}, "cv2-magsac": {"px_th": 4.0, "ratio_th": null, "source": "run.py tune h (val)"}, "poselib": {"px_th": 8.0, "ratio_th": null, "source": "run.py tune h (val)"}, "poselib-prosac": {"px_th": 4.0, "ratio_th": null, "source": "run.py tune h (val), after score-orientation fix"}, "pydegensac": {"px_th": 1.0, "ratio_th": null, "source": "run.py tune h (val)"}}, "HPatchesSeq": {"cv2-ransac": {"px_th": 16.0, "ratio_th": null, "source": "run.py tune h (val)"}, "cv2-magsac": {"px_th": 64.0, "ratio_th": null, "source": "run.py tune h (val)"}, "poselib": {"px_th": 16.0, "ratio_th": null, "source": "run.py tune h (val)"}, "poselib-prosac": {"px_th": 16.0, "ratio_th": null, "source": "run.py tune h (val), after score-orientation fix"}, "pydegensac": {"px_th": 4.0, "ratio_th": null, "source": "run.py tune h (val)"}}}`
- `base:m1-master+lapackfix` — h_m1_basefix2.jsonl, pydegensac 0.3.0, cv2 4.13.0, poselib 2.0.5, repo aa111c8, host Dmytros-MacBook-Air.local
  - config: `{"EVD": {"pydegensac": {"px_th": 1.0, "ratio_th": null, "source": "run.py tune h (val)"}}, "HPatchesSeq": {"pydegensac": {"px_th": 4.0, "ratio_th": null, "source": "run.py tune h (val)"}}}`

### EVD — best per method

| method | best mAA | 95% CI | at budget | mean ms/pair | d mAA vs poselib-prosac (paired) |
|---|---|---|---|---|---|
| poselib-prosac | 0.4500 | 0.2000-0.6875 | 100 | 0.84 | leader |
| poselib | 0.4375 | 0.1875-0.6875 | 25000 | 2.95 | -0.0126 (-0.0375, +0.0000) |
| cv2-magsac | 0.4125 | 0.1750-0.6500 | 6400 | 1.05 | -0.0379 (-0.0750, -0.0125) * |
| pydegensac (branch) | 0.4125 | 0.1750-0.6000 | 25000 | 7.79 | -0.0374 (-0.1125, +0.0375) |
| pydegensac (base) | 0.4000 | 0.1750-0.6375 | 25000 | 31.35 | -0.0500 (-0.0875, -0.0125) * |
| cv2-ransac | 0.3375 | 0.0994-0.6125 | 1600 | 18.76 | -0.1131 (-0.2500, +0.0000) |

`*` = paired CI excludes zero. Each method is taken at its own best budget, which flatters every method equally.

#### pydegensac: base vs branch, same budget

| budget | mAA base | mAA branch | d mAA (95% CI, paired) | ms base | ms branch | speedup |
|---|---|---|---|---|---|---|
| 10 | 0.1000 | 0.0750 | -0.0253 (-0.0750, +0.0000) | 0.35 | 0.25 | 1.39x |
| 25 | 0.2125 | 0.1000 | -0.1135 (-0.2875, +0.0000) | 0.34 | 0.28 | 1.21x |
| 100 | 0.2875 | 0.2125 | -0.0756 (-0.2125, +0.0250) | 0.95 | 0.52 | 1.84x |
| 400 | 0.2875 | 0.2875 | -0.0004 (-0.0750, +0.0750) | 1.78 | 0.92 | 1.93x |
| 1600 | 0.2875 | 0.2625 | -0.0257 (-0.0875, +0.0500) | 4.01 | 1.58 | 2.53x |
| 6400 | 0.3625 | 0.3375 | -0.0253 (-0.0875, +0.0250) | 10.15 | 3.39 | 2.99x |
| 25000 | 0.4000 | 0.4125 | +0.0126 (-0.0750, +0.1000) | 31.35 | 7.79 | 4.02x |
| **all** | | | | 48.91 | 14.74 | **3.32x** |

### HPatchesSeq — best per method

| method | best mAA | 95% CI | at budget | mean ms/pair | d mAA vs poselib-prosac (paired) |
|---|---|---|---|---|---|
| poselib-prosac | 0.9297 | 0.8903-0.9628 | 1600 | 7.97 | leader |
| cv2-magsac | 0.9269 | 0.8896-0.9586 | 25000 | 0.96 | -0.0027 (-0.0124, +0.0076) |
| poselib | 0.9262 | 0.8876-0.9593 | 25000 | 5.22 | -0.0034 (-0.0214, +0.0097) |
| pydegensac (branch) | 0.9200 | 0.8814-0.9538 | 25000 | 2.56 | -0.0095 (-0.0152, -0.0041) * |
| cv2-ransac | 0.9103 | 0.8710-0.9448 | 25000 | 28.24 | -0.0192 (-0.0297, -0.0090) * |
| pydegensac (base) | 0.9041 | 0.8621-0.9393 | 1600 | 4.39 | -0.0253 (-0.0434, -0.0117) * |

`*` = paired CI excludes zero. Each method is taken at its own best budget, which flatters every method equally.

#### pydegensac: base vs branch, same budget

| budget | mAA base | mAA branch | d mAA (95% CI, paired) | ms base | ms branch | speedup |
|---|---|---|---|---|---|---|
| 10 | 0.7048 | 0.7228 | +0.0179 (-0.0524, +0.0848) | 2.79 | 0.98 | 2.85x |
| 25 | 0.7641 | 0.7938 | +0.0304 (-0.0159, +0.0779) | 2.94 | 1.08 | 2.73x |
| 100 | 0.8531 | 0.8366 | -0.0167 (-0.0579, +0.0235) | 3.32 | 1.21 | 2.74x |
| 400 | 0.8862 | 0.8945 | +0.0079 (-0.0110, +0.0269) | 3.78 | 1.43 | 2.63x |
| 1600 | 0.9041 | 0.9069 | +0.0025 (-0.0145, +0.0172) | 4.39 | 1.57 | 2.79x |
| 6400 | 0.9028 | 0.9131 | +0.0103 (+0.0007, +0.0214) * | 5.69 | 1.94 | 2.93x |
| 25000 | 0.9028 | 0.9200 | +0.0171 (+0.0028, +0.0359) * | 8.43 | 2.56 | 3.30x |
| **all** | | | | 31.33 | 10.78 | **2.91x** |
wrote results/time_maa_h_m1.png
