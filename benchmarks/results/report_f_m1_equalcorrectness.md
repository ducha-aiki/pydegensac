# Fundamental time-mAA

- `branch:m1-f349a6c` — f_m1_branch.jsonl, pydegensac 0.3.0, cv2 4.13.0, poselib 2.0.5, repo f349a6c, host Dmytros-MacBook-Air.local
  - config: `{"st_peters_square": {"cv2-ransac": {"px_th": 0.5, "ratio_th": 0.85, "source": "run.py tune f (300-pair tuning subset)"}, "cv2-magsac": {"px_th": 0.25, "ratio_th": 0.85, "source": "run.py tune f (300-pair tuning subset)"}, "poselib": {"px_th": 0.5, "ratio_th": 0.85, "source": "run.py tune f (300-pair tuning subset)"}, "poselib-prosac": {"px_th": 0.5, "ratio_th": 0.85, "source": "run.py tune f; ratio 0.85 not 0.90 (tie-break, see below)"}, "pydegensac": {"px_th": 0.5, "ratio_th": 0.8, "source": "run.py tune f (300-pair tuning subset)"}}}`
- `base:m1-master+lapackfix` — f_m1_basefix.jsonl, pydegensac 0.3.0, cv2 4.13.0, poselib 2.0.5, repo f349a6c, host Dmytros-MacBook-Air.local
  - config: `{"st_peters_square": {"pydegensac": {"px_th": 0.5, "ratio_th": 0.8, "source": "run.py tune f (300-pair tuning subset)"}}}`

### st_peters_square — best per method

| method | best mAA | 95% CI | at budget | mean ms/pair | d mAA vs poselib-prosac (paired) |
|---|---|---|---|---|---|
| poselib-prosac | 0.4570 | 0.4228-0.4873 | 50000 | 50.70 | leader |
| poselib | 0.4392 | 0.4070-0.4698 | 50000 | 49.19 | -0.0178 (-0.0413, +0.0060) |
| pydegensac (base) | 0.4355 | 0.4002-0.4672 | 50000 | 137.31 | -0.0217 (-0.0443, +0.0010) |
| pydegensac (branch) | 0.4218 | 0.3895-0.4533 | 50000 | 96.11 | -0.0351 (-0.0560, -0.0127) * |
| cv2-magsac | 0.3753 | 0.3412-0.4063 | 50000 | 49.65 | -0.0820 (-0.1098, -0.0562) * |
| cv2-ransac | 0.3403 | 0.3103-0.3693 | 50000 | 240.18 | -0.1165 (-0.1425, -0.0898) * |

`*` = paired CI excludes zero. Each method is taken at its own best budget, which flatters every method equally.

#### pydegensac: base vs branch, same budget

| budget | mAA base | mAA branch | d mAA (95% CI, paired) | ms base | ms branch | speedup |
|---|---|---|---|---|---|---|
| 125 | 0.3733 | 0.3573 | -0.0162 (-0.0370, +0.0045) | 6.44 | 6.21 | 1.04x |
| 250 | 0.3890 | 0.3770 | -0.0123 (-0.0323, +0.0067) | 7.53 | 7.08 | 1.06x |
| 500 | 0.4063 | 0.3827 | -0.0241 (-0.0440, -0.0053) * | 9.06 | 8.19 | 1.11x |
| 1000 | 0.4127 | 0.3948 | -0.0181 (-0.0387, +0.0010) | 11.20 | 9.77 | 1.15x |
| 2500 | 0.4110 | 0.4022 | -0.0088 (-0.0275, +0.0107) | 16.34 | 13.52 | 1.21x |
| 5000 | 0.4218 | 0.4140 | -0.0077 (-0.0255, +0.0113) | 24.15 | 18.71 | 1.29x |
| 10000 | 0.4213 | 0.4215 | +0.0008 (-0.0172, +0.0198) | 38.61 | 28.23 | 1.37x |
| 25000 | 0.4273 | 0.4197 | -0.0072 (-0.0265, +0.0123) | 77.30 | 54.81 | 1.41x |
| 50000 | 0.4355 | 0.4218 | -0.0134 (-0.0315, +0.0052) | 137.31 | 96.11 | 1.43x |
| **all** | | | | 327.95 | 242.62 | **1.35x** |
wrote results/time_maa_f_m1_equalcorrectness.png
