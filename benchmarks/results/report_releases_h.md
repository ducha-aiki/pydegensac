# Homography time-mAA

- `local-0.2.2:f7208bc` — rel_h_local-0.2.2.jsonl, pydegensac 0.2.2, cv2 4.13.0, poselib missing, repo e8c028a, host DM-PC
  - config: `{"EVD": {"pydegensac": {"px_th": 1.0, "ratio_th": null, "source": "run.py tune h (val)"}}, "HPatchesSeq": {"pydegensac": {"px_th": 4.0, "ratio_th": null, "source": "run.py tune h (val)"}}}`
- `local-branch:e8c028a` — rel_h_local-branch.jsonl, pydegensac 0.3.0, cv2 4.13.0, poselib missing, repo e8c028a, host DM-PC
  - config: `{"EVD": {"pydegensac": {"px_th": 1.0, "ratio_th": null, "source": "run.py tune h (val)"}}, "HPatchesSeq": {"pydegensac": {"px_th": 4.0, "ratio_th": null, "source": "run.py tune h (val)"}}}`
- `local-master:08464ca` — rel_h_local-master.jsonl, pydegensac 0.3.0, cv2 4.13.0, poselib missing, repo e8c028a, host DM-PC
  - config: `{"EVD": {"pydegensac": {"px_th": 1.0, "ratio_th": null, "source": "run.py tune h (val)"}}, "HPatchesSeq": {"pydegensac": {"px_th": 4.0, "ratio_th": null, "source": "run.py tune h (val)"}}}`
- `pypi-0.1.2:wheel` — rel_h_pypi-0.1.2.jsonl, pydegensac 0.1.2, cv2 4.13.0, poselib missing, repo e8c028a, host DM-PC
  - config: `{"EVD": {"pydegensac": {"px_th": 1.0, "ratio_th": null, "source": "run.py tune h (val)"}}, "HPatchesSeq": {"pydegensac": {"px_th": 4.0, "ratio_th": null, "source": "run.py tune h (val)"}}}`
- `pypi-0.2.1:wheel` — rel_h_pypi-0.2.1.jsonl, pydegensac 0.2.1, cv2 4.13.0, poselib missing, repo e8c028a, host DM-PC
  - config: `{"EVD": {"pydegensac": {"px_th": 1.0, "ratio_th": null, "source": "run.py tune h (val)"}}, "HPatchesSeq": {"pydegensac": {"px_th": 4.0, "ratio_th": null, "source": "run.py tune h (val)"}}}`
- `pypi-0.2.2:wheel` — rel_h_pypi-0.2.2.jsonl, pydegensac 0.2.2, cv2 4.13.0, poselib missing, repo e8c028a, host DM-PC
  - config: `{"EVD": {"pydegensac": {"px_th": 1.0, "ratio_th": null, "source": "run.py tune h (val)"}}, "HPatchesSeq": {"pydegensac": {"px_th": 4.0, "ratio_th": null, "source": "run.py tune h (val)"}}}`

### EVD — best per method

| method | best mAA | 95% CI | at budget | mean ms/pair | d mAA vs pydegensac (local-master) (paired) |
|---|---|---|---|---|---|
| pydegensac (local-master) | 0.4625 | 0.2500-0.6625 | 6400 | 7.78 | leader |
| pydegensac (local-branch) | 0.4125 | 0.1997-0.6375 | 25000 | 13.22 | -0.0501 (-0.2125, +0.1375) |
| pydegensac (pypi-0.1.2) | 0.4125 | 0.1875-0.6500 | 1600 | 5.01 | -0.0504 (-0.1750, +0.0625) |
| pydegensac (pypi-0.2.1) | 0.3875 | 0.1497-0.6375 | 6400 | 10.66 | -0.0758 (-0.1875, +0.0125) |
| pydegensac (pypi-0.2.2) | 0.3875 | 0.1625-0.6253 | 1600 | 5.34 | -0.0758 (-0.2000, +0.0375) |
| pydegensac (local-0.2.2) | 0.3750 | 0.1375-0.6375 | 25000 | 24.43 | -0.0880 (-0.2000, +0.0125) |

`*` = paired CI excludes zero. Each method is taken at its own best budget, which flatters every method equally.


### EVD

| method | budget | mAA | 95% CI | mean ms/pair | median ms/pair | failures |
|---|---|---|---|---|---|---|
| pydegensac (local-0.2.2) | 10 | 0.0750 | 0.0000-0.2250 | 0.29 | 0.25 | 5/8 |
| pydegensac (local-0.2.2) | 25 | 0.0750 | 0.0000-0.2000 | 0.33 | 0.22 | 2/8 |
| pydegensac (local-0.2.2) | 100 | 0.2625 | 0.0625-0.5000 | 0.69 | 0.29 | 0/8 |
| pydegensac (local-0.2.2) | 400 | 0.2875 | 0.0625-0.5750 | 1.61 | 0.79 | 0/8 |
| pydegensac (local-0.2.2) | 1600 | 0.3625 | 0.1250-0.6250 | 3.47 | 3.54 | 0/8 |
| pydegensac (local-0.2.2) | 6400 | 0.3625 | 0.1250-0.6250 | 8.66 | 8.51 | 0/8 |
| pydegensac (local-0.2.2) | 25000 | 0.3750 | 0.1375-0.6375 | 24.43 | 26.98 | 0/8 |
| pydegensac (local-branch) | 10 | 0.1000 | 0.0000-0.3000 | 0.28 | 0.15 | 2/8 |
| pydegensac (local-branch) | 25 | 0.1125 | 0.0000-0.3375 | 0.30 | 0.15 | 1/8 |
| pydegensac (local-branch) | 100 | 0.2250 | 0.0375-0.4625 | 0.72 | 0.37 | 0/8 |
| pydegensac (local-branch) | 400 | 0.2375 | 0.0375-0.5000 | 1.19 | 0.86 | 0/8 |
| pydegensac (local-branch) | 1600 | 0.2625 | 0.0625-0.5000 | 2.42 | 1.40 | 0/8 |
| pydegensac (local-branch) | 6400 | 0.3375 | 0.1000-0.5750 | 5.28 | 4.45 | 0/8 |
| pydegensac (local-branch) | 25000 | 0.4125 | 0.1997-0.6375 | 13.22 | 10.24 | 0/8 |
| pydegensac (local-master) | 10 | 0.0875 | 0.0000-0.2625 | 0.56 | 0.22 | 4/8 |
| pydegensac (local-master) | 25 | 0.1000 | 0.0000-0.3000 | 0.29 | 0.15 | 2/8 |
| pydegensac (local-master) | 100 | 0.1625 | 0.0000-0.3375 | 0.53 | 0.27 | 1/8 |
| pydegensac (local-master) | 400 | 0.2875 | 0.0750-0.5500 | 1.60 | 1.21 | 0/8 |
| pydegensac (local-master) | 1600 | 0.4500 | 0.2125-0.6750 | 3.51 | 3.01 | 0/8 |
| pydegensac (local-master) | 6400 | 0.4625 | 0.2500-0.6625 | 7.78 | 7.41 | 0/8 |
| pydegensac (local-master) | 25000 | 0.4625 | 0.2500-0.6625 | 23.32 | 24.51 | 0/8 |
| pydegensac (pypi-0.1.2) | 10 | 0.1125 | 0.0000-0.3375 | 0.28 | 0.19 | 3/8 |
| pydegensac (pypi-0.1.2) | 25 | 0.0750 | 0.0000-0.2250 | 0.26 | 0.17 | 0/8 |
| pydegensac (pypi-0.1.2) | 100 | 0.1750 | 0.0000-0.4000 | 0.94 | 0.77 | 0/8 |
| pydegensac (pypi-0.1.2) | 400 | 0.3375 | 0.1116-0.6000 | 2.50 | 1.89 | 0/8 |
| pydegensac (pypi-0.1.2) | 1600 | 0.4125 | 0.1875-0.6500 | 5.01 | 4.18 | 0/8 |
| pydegensac (pypi-0.1.2) | 6400 | 0.4125 | 0.1875-0.6500 | 10.85 | 9.46 | 0/8 |
| pydegensac (pypi-0.1.2) | 25000 | 0.4125 | 0.1875-0.6500 | 27.33 | 25.31 | 0/8 |
| pydegensac (pypi-0.2.1) | 10 | 0.0750 | 0.0000-0.2250 | 0.30 | 0.17 | 1/8 |
| pydegensac (pypi-0.2.1) | 25 | 0.1625 | 0.0000-0.3875 | 0.37 | 0.24 | 0/8 |
| pydegensac (pypi-0.2.1) | 100 | 0.2375 | 0.0375-0.5000 | 0.85 | 0.78 | 0/8 |
| pydegensac (pypi-0.2.1) | 400 | 0.2250 | 0.0250-0.5000 | 1.77 | 1.41 | 0/8 |
| pydegensac (pypi-0.2.1) | 1600 | 0.2875 | 0.0875-0.5625 | 4.24 | 3.37 | 0/8 |
| pydegensac (pypi-0.2.1) | 6400 | 0.3875 | 0.1497-0.6375 | 10.66 | 7.89 | 0/8 |
| pydegensac (pypi-0.2.1) | 25000 | 0.3875 | 0.1497-0.6375 | 26.76 | 22.78 | 0/8 |
| pydegensac (pypi-0.2.2) | 10 | 0.1875 | 0.0000-0.4500 | 0.36 | 0.27 | 5/8 |
| pydegensac (pypi-0.2.2) | 25 | 0.1750 | 0.0000-0.4000 | 0.31 | 0.15 | 5/8 |
| pydegensac (pypi-0.2.2) | 100 | 0.2375 | 0.0500-0.4750 | 0.82 | 0.48 | 0/8 |
| pydegensac (pypi-0.2.2) | 400 | 0.2625 | 0.0500-0.4753 | 2.16 | 1.59 | 0/8 |
| pydegensac (pypi-0.2.2) | 1600 | 0.3875 | 0.1625-0.6253 | 5.34 | 5.18 | 0/8 |
| pydegensac (pypi-0.2.2) | 6400 | 0.3625 | 0.1372-0.6128 | 12.19 | 13.26 | 0/8 |
| pydegensac (pypi-0.2.2) | 25000 | 0.3625 | 0.1372-0.6128 | 33.07 | 34.44 | 0/8 |

### HPatchesSeq — best per method

| method | best mAA | 95% CI | at budget | mean ms/pair | d mAA vs pydegensac (pypi-0.2.1) (paired) |
|---|---|---|---|---|---|
| pydegensac (pypi-0.2.1) | 0.9193 | 0.8800-0.9538 | 25000 | 7.43 | leader |
| pydegensac (local-branch) | 0.9172 | 0.8779-0.9510 | 25000 | 4.43 | -0.0020 (-0.0103, +0.0062) |
| pydegensac (local-0.2.2) | 0.9138 | 0.8752-0.9476 | 25000 | 6.24 | -0.0054 (-0.0152, +0.0041) |
| pydegensac (pypi-0.2.2) | 0.9124 | 0.8738-0.9469 | 6400 | 5.36 | -0.0070 (-0.0159, +0.0007) |
| pydegensac (pypi-0.1.2) | 0.9103 | 0.8697-0.9448 | 25000 | 7.69 | -0.0090 (-0.0262, +0.0048) |
| pydegensac (local-master) | 0.9034 | 0.8628-0.9386 | 6400 | 4.11 | -0.0159 (-0.0338, -0.0014) * |

`*` = paired CI excludes zero. Each method is taken at its own best budget, which flatters every method equally.


### HPatchesSeq

| method | budget | mAA | 95% CI | mean ms/pair | median ms/pair | failures |
|---|---|---|---|---|---|---|
| pydegensac (local-0.2.2) | 10 | 0.7014 | 0.6255-0.7703 | 1.56 | 1.77 | 21/145 |
| pydegensac (local-0.2.2) | 25 | 0.7841 | 0.7186-0.8435 | 1.74 | 1.87 | 6/145 |
| pydegensac (local-0.2.2) | 100 | 0.8579 | 0.8055-0.9055 | 2.08 | 2.14 | 0/145 |
| pydegensac (local-0.2.2) | 400 | 0.8759 | 0.8276-0.9186 | 2.43 | 2.29 | 0/145 |
| pydegensac (local-0.2.2) | 1600 | 0.8924 | 0.8483-0.9310 | 2.84 | 2.43 | 0/145 |
| pydegensac (local-0.2.2) | 6400 | 0.9131 | 0.8745-0.9469 | 3.98 | 2.51 | 0/145 |
| pydegensac (local-0.2.2) | 25000 | 0.9138 | 0.8752-0.9476 | 6.24 | 2.54 | 0/145 |
| pydegensac (local-branch) | 10 | 0.7345 | 0.6627-0.8007 | 1.58 | 1.79 | 24/145 |
| pydegensac (local-branch) | 25 | 0.7848 | 0.7207-0.8455 | 1.71 | 1.89 | 5/145 |
| pydegensac (local-branch) | 100 | 0.8400 | 0.7841-0.8938 | 1.98 | 2.07 | 0/145 |
| pydegensac (local-branch) | 400 | 0.8745 | 0.8241-0.9207 | 2.35 | 2.24 | 0/145 |
| pydegensac (local-branch) | 1600 | 0.8993 | 0.8552-0.9400 | 2.69 | 2.38 | 0/145 |
| pydegensac (local-branch) | 6400 | 0.9048 | 0.8614-0.9421 | 3.35 | 2.51 | 0/145 |
| pydegensac (local-branch) | 25000 | 0.9172 | 0.8779-0.9510 | 4.43 | 2.52 | 0/145 |
| pydegensac (local-master) | 10 | 0.7179 | 0.6421-0.7855 | 1.63 | 1.78 | 18/145 |
| pydegensac (local-master) | 25 | 0.7959 | 0.7317-0.8545 | 1.78 | 1.85 | 5/145 |
| pydegensac (local-master) | 100 | 0.8545 | 0.8014-0.9028 | 2.06 | 2.08 | 0/145 |
| pydegensac (local-master) | 400 | 0.8738 | 0.8248-0.9172 | 2.54 | 2.31 | 0/145 |
| pydegensac (local-master) | 1600 | 0.8848 | 0.8386-0.9248 | 3.01 | 2.42 | 0/145 |
| pydegensac (local-master) | 6400 | 0.9034 | 0.8628-0.9386 | 4.11 | 2.53 | 0/145 |
| pydegensac (local-master) | 25000 | 0.9034 | 0.8628-0.9386 | 6.34 | 2.46 | 0/145 |
| pydegensac (pypi-0.1.2) | 10 | 0.7090 | 0.6317-0.7834 | 2.00 | 2.29 | 27/145 |
| pydegensac (pypi-0.1.2) | 25 | 0.7800 | 0.7110-0.8407 | 2.22 | 2.46 | 11/145 |
| pydegensac (pypi-0.1.2) | 100 | 0.8290 | 0.7683-0.8855 | 2.62 | 2.71 | 0/145 |
| pydegensac (pypi-0.1.2) | 400 | 0.8800 | 0.8324-0.9241 | 3.22 | 2.92 | 0/145 |
| pydegensac (pypi-0.1.2) | 1600 | 0.9021 | 0.8586-0.9400 | 3.84 | 3.09 | 0/145 |
| pydegensac (pypi-0.1.2) | 6400 | 0.9055 | 0.8627-0.9434 | 5.24 | 3.35 | 0/145 |
| pydegensac (pypi-0.1.2) | 25000 | 0.9103 | 0.8697-0.9448 | 7.69 | 3.38 | 0/145 |
| pydegensac (pypi-0.2.1) | 10 | 0.7434 | 0.6703-0.8104 | 2.06 | 2.33 | 17/145 |
| pydegensac (pypi-0.2.1) | 25 | 0.7986 | 0.7338-0.8593 | 2.29 | 2.53 | 7/145 |
| pydegensac (pypi-0.2.1) | 100 | 0.8710 | 0.8186-0.9172 | 2.60 | 2.71 | 0/145 |
| pydegensac (pypi-0.2.1) | 400 | 0.8848 | 0.8338-0.9262 | 3.18 | 2.96 | 0/145 |
| pydegensac (pypi-0.2.1) | 1600 | 0.9000 | 0.8552-0.9379 | 3.77 | 3.05 | 0/145 |
| pydegensac (pypi-0.2.1) | 6400 | 0.9179 | 0.8793-0.9524 | 5.08 | 3.25 | 0/145 |
| pydegensac (pypi-0.2.1) | 25000 | 0.9193 | 0.8800-0.9538 | 7.43 | 3.25 | 0/145 |
| pydegensac (pypi-0.2.2) | 10 | 0.7152 | 0.6414-0.7821 | 2.06 | 2.36 | 16/145 |
| pydegensac (pypi-0.2.2) | 25 | 0.7972 | 0.7331-0.8559 | 2.28 | 2.52 | 3/145 |
| pydegensac (pypi-0.2.2) | 100 | 0.8593 | 0.8062-0.9076 | 2.72 | 2.77 | 0/145 |
| pydegensac (pypi-0.2.2) | 400 | 0.8759 | 0.8269-0.9200 | 3.27 | 2.95 | 0/145 |
| pydegensac (pypi-0.2.2) | 1600 | 0.8890 | 0.8414-0.9310 | 3.93 | 3.13 | 0/145 |
| pydegensac (pypi-0.2.2) | 6400 | 0.9124 | 0.8738-0.9469 | 5.36 | 3.24 | 0/145 |
| pydegensac (pypi-0.2.2) | 25000 | 0.9117 | 0.8724-0.9462 | 8.29 | 3.32 | 0/145 |
wrote results/releases_h.png
