# Fundamental time-mAA

- `local-0.2.2:f7208bc` — rel_f_local-0.2.2.jsonl, pydegensac 0.2.2, cv2 4.13.0, poselib missing, repo e8c028a, host DM-PC
  - config: `{"st_peters_square": {"pydegensac": {"px_th": 0.5, "ratio_th": 0.8, "source": "run.py tune f (300-pair tuning subset)"}}}`
- `local-branch:e8c028a` — rel_f_local-branch.jsonl, pydegensac 0.3.0, cv2 4.13.0, poselib missing, repo e8c028a, host DM-PC
  - config: `{"st_peters_square": {"pydegensac": {"px_th": 0.5, "ratio_th": 0.8, "source": "run.py tune f (300-pair tuning subset)"}}}`
- `local-master:08464ca` — rel_f_local-master.jsonl, pydegensac 0.3.0, cv2 4.13.0, poselib missing, repo e8c028a, host DM-PC
  - config: `{"st_peters_square": {"pydegensac": {"px_th": 0.5, "ratio_th": 0.8, "source": "run.py tune f (300-pair tuning subset)"}}}`
- `pypi-0.1.2:wheel` — rel_f_pypi-0.1.2.jsonl, pydegensac 0.1.2, cv2 4.13.0, poselib missing, repo e8c028a, host DM-PC
  - config: `{"st_peters_square": {"pydegensac": {"px_th": 0.5, "ratio_th": 0.8, "source": "run.py tune f (300-pair tuning subset)"}}}`
- `pypi-0.2.1:wheel` — rel_f_pypi-0.2.1.jsonl, pydegensac 0.2.1, cv2 4.13.0, poselib missing, repo e8c028a, host DM-PC
  - config: `{"st_peters_square": {"pydegensac": {"px_th": 0.5, "ratio_th": 0.8, "source": "run.py tune f (300-pair tuning subset)"}}}`
- `pypi-0.2.2:wheel` — rel_f_pypi-0.2.2.jsonl, pydegensac 0.2.2, cv2 4.13.0, poselib missing, repo e8c028a, host DM-PC
  - config: `{"st_peters_square": {"pydegensac": {"px_th": 0.5, "ratio_th": 0.8, "source": "run.py tune f (300-pair tuning subset)"}}}`

### st_peters_square — best per method

| method | best mAA | 95% CI | at budget | mean ms/pair | d mAA vs pydegensac (pypi-0.2.1) (paired) |
|---|---|---|---|---|---|
| pydegensac (pypi-0.2.1) | 0.4390 | 0.4042-0.4697 | 50000 | 131.66 | leader |
| pydegensac (pypi-0.2.2) | 0.4350 | 0.4012-0.4658 | 25000 | 78.92 | -0.0038 (-0.0233, +0.0140) |
| pydegensac (local-master) | 0.4320 | 0.3980-0.4640 | 50000 | 134.52 | -0.0068 (-0.0262, +0.0118) |
| pydegensac (local-branch) | 0.4285 | 0.3943-0.4612 | 10000 | 30.72 | -0.0104 (-0.0313, +0.0093) |
| pydegensac (pypi-0.1.2) | 0.4267 | 0.3927-0.4590 | 50000 | 129.74 | -0.0122 (-0.0320, +0.0077) |
| pydegensac (local-0.2.2) | 0.4233 | 0.3898-0.4548 | 50000 | 131.44 | -0.0154 (-0.0355, +0.0047) |

`*` = paired CI excludes zero. Each method is taken at its own best budget, which flatters every method equally.


### st_peters_square

| method | budget | mAA | 95% CI | mean ms/pair | median ms/pair | failures |
|---|---|---|---|---|---|---|
| pydegensac (local-0.2.2) | 125 | 0.3743 | 0.3407-0.4053 | 5.55 | 4.42 | 0/600 |
| pydegensac (local-0.2.2) | 250 | 0.3787 | 0.3447-0.4107 | 6.54 | 5.23 | 0/600 |
| pydegensac (local-0.2.2) | 500 | 0.3897 | 0.3552-0.4213 | 7.83 | 6.16 | 0/600 |
| pydegensac (local-0.2.2) | 1000 | 0.4023 | 0.3692-0.4335 | 9.86 | 7.87 | 0/600 |
| pydegensac (local-0.2.2) | 2500 | 0.4153 | 0.3815-0.4463 | 14.79 | 11.28 | 0/600 |
| pydegensac (local-0.2.2) | 5000 | 0.4225 | 0.3873-0.4545 | 22.04 | 16.34 | 0/600 |
| pydegensac (local-0.2.2) | 10000 | 0.4217 | 0.3878-0.4538 | 35.53 | 25.11 | 0/600 |
| pydegensac (local-0.2.2) | 25000 | 0.4208 | 0.3867-0.4525 | 72.94 | 48.73 | 0/600 |
| pydegensac (local-0.2.2) | 50000 | 0.4233 | 0.3898-0.4548 | 131.44 | 86.16 | 0/600 |
| pydegensac (local-branch) | 125 | 0.3708 | 0.3372-0.4033 | 5.40 | 4.06 | 0/600 |
| pydegensac (local-branch) | 250 | 0.3850 | 0.3518-0.4173 | 6.29 | 4.84 | 0/600 |
| pydegensac (local-branch) | 500 | 0.3897 | 0.3553-0.4215 | 7.50 | 5.80 | 0/600 |
| pydegensac (local-branch) | 1000 | 0.3942 | 0.3595-0.4252 | 9.32 | 7.09 | 0/600 |
| pydegensac (local-branch) | 2500 | 0.3993 | 0.3643-0.4303 | 13.38 | 9.43 | 0/600 |
| pydegensac (local-branch) | 5000 | 0.4208 | 0.3867-0.4522 | 19.32 | 12.98 | 0/600 |
| pydegensac (local-branch) | 10000 | 0.4285 | 0.3943-0.4612 | 30.72 | 19.25 | 0/600 |
| pydegensac (local-branch) | 25000 | 0.4258 | 0.3902-0.4580 | 62.19 | 35.74 | 0/600 |
| pydegensac (local-branch) | 50000 | 0.4203 | 0.3860-0.4517 | 111.28 | 60.22 | 0/600 |
| pydegensac (local-master) | 125 | 0.3633 | 0.3297-0.3950 | 5.60 | 4.35 | 0/600 |
| pydegensac (local-master) | 250 | 0.3758 | 0.3422-0.4058 | 6.53 | 4.92 | 0/600 |
| pydegensac (local-master) | 500 | 0.3830 | 0.3498-0.4137 | 7.83 | 6.21 | 0/600 |
| pydegensac (local-master) | 1000 | 0.3923 | 0.3592-0.4228 | 9.89 | 7.75 | 0/600 |
| pydegensac (local-master) | 2500 | 0.4035 | 0.3702-0.4340 | 14.86 | 11.06 | 0/600 |
| pydegensac (local-master) | 5000 | 0.4128 | 0.3800-0.4427 | 22.20 | 15.90 | 0/600 |
| pydegensac (local-master) | 10000 | 0.4185 | 0.3853-0.4497 | 35.82 | 24.63 | 0/600 |
| pydegensac (local-master) | 25000 | 0.4208 | 0.3867-0.4532 | 73.92 | 47.55 | 0/600 |
| pydegensac (local-master) | 50000 | 0.4320 | 0.3980-0.4640 | 134.52 | 87.80 | 0/600 |
| pydegensac (pypi-0.1.2) | 125 | 0.3768 | 0.3447-0.4093 | 6.98 | 5.20 | 0/600 |
| pydegensac (pypi-0.1.2) | 250 | 0.3870 | 0.3547-0.4195 | 8.05 | 6.02 | 0/600 |
| pydegensac (pypi-0.1.2) | 500 | 0.3933 | 0.3588-0.4252 | 9.43 | 7.24 | 0/600 |
| pydegensac (pypi-0.1.2) | 1000 | 0.4032 | 0.3692-0.4348 | 11.57 | 8.89 | 0/600 |
| pydegensac (pypi-0.1.2) | 2500 | 0.4130 | 0.3797-0.4447 | 16.52 | 12.10 | 0/600 |
| pydegensac (pypi-0.1.2) | 5000 | 0.4195 | 0.3862-0.4507 | 23.57 | 16.36 | 0/600 |
| pydegensac (pypi-0.1.2) | 10000 | 0.4230 | 0.3890-0.4542 | 36.51 | 24.22 | 0/600 |
| pydegensac (pypi-0.1.2) | 25000 | 0.4208 | 0.3870-0.4528 | 73.13 | 47.03 | 0/600 |
| pydegensac (pypi-0.1.2) | 50000 | 0.4267 | 0.3927-0.4590 | 129.74 | 82.16 | 0/600 |
| pydegensac (pypi-0.2.1) | 125 | 0.3690 | 0.3350-0.4000 | 7.16 | 5.66 | 0/600 |
| pydegensac (pypi-0.2.1) | 250 | 0.3857 | 0.3527-0.4160 | 8.34 | 6.70 | 0/600 |
| pydegensac (pypi-0.2.1) | 500 | 0.3953 | 0.3628-0.4255 | 9.72 | 7.74 | 0/600 |
| pydegensac (pypi-0.2.1) | 1000 | 0.4078 | 0.3740-0.4383 | 11.86 | 9.63 | 0/600 |
| pydegensac (pypi-0.2.1) | 2500 | 0.4187 | 0.3848-0.4502 | 16.77 | 12.71 | 0/600 |
| pydegensac (pypi-0.2.1) | 5000 | 0.4240 | 0.3897-0.4567 | 23.94 | 17.19 | 0/600 |
| pydegensac (pypi-0.2.1) | 10000 | 0.4300 | 0.3952-0.4618 | 36.97 | 25.18 | 0/600 |
| pydegensac (pypi-0.2.1) | 25000 | 0.4292 | 0.3937-0.4593 | 73.50 | 46.73 | 0/600 |
| pydegensac (pypi-0.2.1) | 50000 | 0.4390 | 0.4042-0.4697 | 131.66 | 81.43 | 0/600 |
| pydegensac (pypi-0.2.2) | 125 | 0.3783 | 0.3440-0.4120 | 7.08 | 5.64 | 0/600 |
| pydegensac (pypi-0.2.2) | 250 | 0.3933 | 0.3603-0.4263 | 8.26 | 6.68 | 0/600 |
| pydegensac (pypi-0.2.2) | 500 | 0.4000 | 0.3672-0.4337 | 9.79 | 7.81 | 0/600 |
| pydegensac (pypi-0.2.2) | 1000 | 0.4137 | 0.3802-0.4453 | 11.98 | 9.41 | 0/600 |
| pydegensac (pypi-0.2.2) | 2500 | 0.4175 | 0.3850-0.4502 | 17.23 | 13.34 | 0/600 |
| pydegensac (pypi-0.2.2) | 5000 | 0.4163 | 0.3843-0.4473 | 24.79 | 18.32 | 0/600 |
| pydegensac (pypi-0.2.2) | 10000 | 0.4222 | 0.3885-0.4528 | 39.00 | 28.21 | 0/600 |
| pydegensac (pypi-0.2.2) | 25000 | 0.4350 | 0.4012-0.4658 | 78.92 | 52.06 | 0/600 |
| pydegensac (pypi-0.2.2) | 50000 | 0.4292 | 0.3958-0.4613 | 141.21 | 94.01 | 0/600 |
wrote results/releases_f.png
