# Fundamental time-mAA

- `branch:m1-opt@15445ef` — f_m1_opt.jsonl, pydegensac 0.3.0, cv2 4.13.0, poselib 2.0.5, repo 15445ef, host Dmytros-MacBook-Air.local
  - config: `{"st_peters_square": {"cv2-ransac": {"px_th": 0.5, "ratio_th": 0.85, "source": "run.py tune f (300-pair tuning subset)"}, "cv2-magsac": {"px_th": 0.25, "ratio_th": 0.85, "source": "run.py tune f (300-pair tuning subset)"}, "poselib": {"px_th": 0.5, "ratio_th": 0.85, "source": "run.py tune f (300-pair tuning subset)"}, "poselib-prosac": {"px_th": 0.5, "ratio_th": 0.85, "source": "run.py tune f; ratio 0.85 not 0.90 (tie-break, see below)"}, "pydegensac": {"px_th": 0.5, "ratio_th": 0.8, "source": "run.py tune f (300-pair tuning subset)"}}}`
- `base:m1-master+lapackfix` — f_m1_basefix2.jsonl, pydegensac 0.3.0, cv2 4.13.0, poselib 2.0.5, repo aa111c8, host Dmytros-MacBook-Air.local
  - config: `{"st_peters_square": {"pydegensac": {"px_th": 0.5, "ratio_th": 0.8, "source": "run.py tune f (300-pair tuning subset)"}}}`

### st_peters_square — best per method

| method | best mAA | 95% CI | at budget | mean ms/pair | d mAA vs poselib-prosac (paired) |
|---|---|---|---|---|---|
| poselib-prosac | 0.4570 | 0.4228-0.4873 | 50000 | 51.81 | leader |
| poselib | 0.4392 | 0.4070-0.4698 | 50000 | 49.75 | -0.0178 (-0.0413, +0.0060) |
| pydegensac (branch) | 0.4358 | 0.4017-0.4670 | 50000 | 49.96 | -0.0211 (-0.0442, +0.0018) |
| pydegensac (base) | 0.4245 | 0.3912-0.4563 | 10000 | 38.75 | -0.0325 (-0.0548, -0.0093) * |
| cv2-magsac | 0.3753 | 0.3412-0.4063 | 50000 | 50.65 | -0.0820 (-0.1098, -0.0562) * |
| cv2-ransac | 0.3403 | 0.3103-0.3693 | 50000 | 243.66 | -0.1165 (-0.1425, -0.0898) * |

`*` = paired CI excludes zero. Each method is taken at its own best budget, which flatters every method equally.

#### pydegensac: base vs branch, same budget

| budget | mAA base | mAA branch | d mAA (95% CI, paired) | ms base | ms branch | speedup |
|---|---|---|---|---|---|---|
| 125 | 0.3802 | 0.3767 | -0.0032 (-0.0243, +0.0175) | 6.30 | 5.03 | 1.25x |
| 250 | 0.3895 | 0.3908 | +0.0019 (-0.0193, +0.0220) | 7.47 | 5.66 | 1.32x |
| 500 | 0.4002 | 0.4002 | +0.0003 (-0.0213, +0.0203) | 8.92 | 6.41 | 1.39x |
| 1000 | 0.4100 | 0.4090 | -0.0004 (-0.0218, +0.0205) | 11.30 | 7.54 | 1.50x |
| 2500 | 0.4130 | 0.4282 | +0.0155 (-0.0052, +0.0355) | 16.49 | 9.65 | 1.71x |
| 5000 | 0.4202 | 0.4287 | +0.0088 (-0.0103, +0.0285) | 24.43 | 12.46 | 1.96x |
| 10000 | 0.4245 | 0.4343 | +0.0099 (-0.0090, +0.0292) | 38.75 | 17.45 | 2.22x |
| 25000 | 0.4240 | 0.4265 | +0.0025 (-0.0163, +0.0218) | 78.29 | 30.60 | 2.56x |
| 50000 | 0.4243 | 0.4358 | +0.0115 (-0.0070, +0.0310) | 138.89 | 49.96 | 2.78x |
| **all** | | | | 330.83 | 144.77 | **2.29x** |
wrote results/time_maa_f_m1.png
