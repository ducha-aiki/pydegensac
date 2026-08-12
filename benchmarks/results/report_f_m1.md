# Fundamental time-mAA

- `branch:m1-f349a6c` — f_m1_branch.jsonl, pydegensac 0.3.0, cv2 4.13.0, poselib 2.0.5, repo f349a6c, host Dmytros-MacBook-Air.local
  - config: `{"st_peters_square": {"cv2-ransac": {"px_th": 0.5, "ratio_th": 0.85, "source": "run.py tune f (300-pair tuning subset)"}, "cv2-magsac": {"px_th": 0.25, "ratio_th": 0.85, "source": "run.py tune f (300-pair tuning subset)"}, "poselib": {"px_th": 0.5, "ratio_th": 0.85, "source": "run.py tune f (300-pair tuning subset)"}, "poselib-prosac": {"px_th": 0.5, "ratio_th": 0.85, "source": "run.py tune f; ratio 0.85 not 0.90 (tie-break, see below)"}, "pydegensac": {"px_th": 0.5, "ratio_th": 0.8, "source": "run.py tune f (300-pair tuning subset)"}}}`
- `base:m1-master` — f_m1_base.jsonl, pydegensac 0.3.0, cv2 4.13.0, poselib 2.0.5, repo f349a6c, host Dmytros-MacBook-Air.local
  - config: `{"st_peters_square": {"pydegensac": {"px_th": 0.5, "ratio_th": 0.8, "source": "run.py tune f (300-pair tuning subset)"}}}`

### st_peters_square — best per method

| method | best mAA | 95% CI | at budget | mean ms/pair | d mAA vs poselib-prosac (paired) |
|---|---|---|---|---|---|
| poselib-prosac | 0.4570 | 0.4228-0.4873 | 50000 | 50.70 | leader |
| poselib | 0.4392 | 0.4070-0.4698 | 50000 | 49.19 | -0.0178 (-0.0413, +0.0060) |
| pydegensac (branch) | 0.4218 | 0.3895-0.4533 | 50000 | 96.11 | -0.0351 (-0.0560, -0.0127) * |
| cv2-magsac | 0.3753 | 0.3412-0.4063 | 50000 | 49.65 | -0.0820 (-0.1098, -0.0562) * |
| cv2-ransac | 0.3403 | 0.3103-0.3693 | 50000 | 240.18 | -0.1165 (-0.1425, -0.0898) * |
| pydegensac (base) | 0.3133 | 0.2803-0.3417 | 50000 | 8.70 | -0.1437 (-0.1723, -0.1142) * |

`*` = paired CI excludes zero. Each method is taken at its own best budget, which flatters every method equally.

#### pydegensac: base vs branch, same budget

| budget | mAA base | mAA branch | d mAA (95% CI, paired) | ms base | ms branch | speedup |
|---|---|---|---|---|---|---|
| 125 | 0.2177 | 0.3573 | +0.1393 (+0.1073, +0.1688) * | 0.45 | 6.21 | 0.07x |
| 250 | 0.2443 | 0.3770 | +0.1323 (+0.1023, +0.1617) * | 0.73 | 7.08 | 0.10x |
| 500 | 0.2642 | 0.3827 | +0.1182 (+0.0905, +0.1457) * | 1.11 | 8.19 | 0.14x |
| 1000 | 0.2807 | 0.3948 | +0.1138 (+0.0877, +0.1417) * | 1.62 | 9.77 | 0.17x |
| 2500 | 0.2977 | 0.4022 | +0.1044 (+0.0778, +0.1317) * | 2.62 | 13.52 | 0.19x |
| 5000 | 0.3058 | 0.4140 | +0.1080 (+0.0820, +0.1350) * | 3.72 | 18.71 | 0.20x |
| 10000 | 0.3092 | 0.4215 | +0.1126 (+0.0855, +0.1397) * | 5.18 | 28.23 | 0.18x |
| 25000 | 0.3120 | 0.4197 | +0.1078 (+0.0805, +0.1347) * | 7.28 | 54.81 | 0.13x |
| 50000 | 0.3133 | 0.4218 | +0.1086 (+0.0813, +0.1365) * | 8.70 | 96.11 | 0.09x |
| **all** | | | | 31.40 | 242.62 | **0.13x** |
wrote results/time_maa_f_m1.png
