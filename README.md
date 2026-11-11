# pydegensac

This repository contains an Python wrapper of RANSAC for homography and fundamental matrix estimation
from sparse correspondences. It implements [LO-RANSAC](https://link.springer.com/chapter/10.1007/978-3-540-45243-0_31) and [DEGENSAC](https://doi.org/10.1109/CVPR.2005.354).

It was originally located in [https://github.com/ducha-aiki/pyransac](https://github.com/ducha-aiki/pyransac), but was renamed to avoid conflict with already existing [pyransac](https://pypi.org/project/pyransac/) in pypi from other author.

> ### macOS users: upgrade
>
> Every macOS build published before **0.3.0** silently skipped all of its
> LAPACK calls. The C sources guarded them behind `#ifdef _WIN32` / `#ifdef
> __linux__`, and macOS defines neither, so the preprocessor removed them: the
> least-squares refits that local optimisation depends on returned an identity
> matrix (F) or an untouched covariance matrix (H). LO-RANSAC's local
> optimisation — the thing that makes this estimator worth using — never ran.
>
> Results were degraded, not broken, so nothing failed loudly: measured on
> public data it cost **0.11 mAA on fundamental and 0.21 on homography**.
> Anything you benchmarked on macOS against other estimators was measuring a
> crippled build. Linux and Windows were unaffected.
>
> To check an existing install:
>
> ```bash
> nm -u $(python -c "import pydegensac,glob,os;print(glob.glob(os.path.dirname(pydegensac.__file__)+'/*.so')[0])") | grep -E 'dgesvd|dsyev'
> ```
>
> Empty output means you have an affected build.

# Performance

Vanilla pydegensac implementation is marginally better than OpenCV one and with degeneracy-check enabled (DEGENSAC) it is the state of the art,
according to the recent study Yin et.al."[Image Matching across Wide Baselines: From Paper to Practice](https://arxiv.org/abs/2003.01587.pdf)", 2020.

![IMW-benchmark](img/ransacs.png)


![IMW-Challenge](img/ransacs2.png)


For homography, pydegensac is worse than newest OpenCV MAGSAC++ (`cv2.USAC_MAGSAC`), but better than OpenCV vanilla RANSAC, according to recent Barath et al. [A Large Scale Homography Benchmark](https://arxiv.org/abs/2302.09997), CVPR2023


![H-benchmark](https://github.com/danini/homography-benchmark/raw/main/assets/heb_benchmark_traditional.png)

## Speed/accuracy against the current field (2026)

`benchmarks/` is a self-contained accuracy-vs-compute benchmark on public data:
fundamental matrix on IMC-2020 PhotoTourism val (600 pairs of
`st_peters_square`, pose mAA 1-10°), homography on EVD (8 test pairs) +
HPatchesSeq (145 test pairs) from the CVPR-2020 RANSAC tutorial data
(reprojection mAA 1-20 px). Every method runs at its own tuned thresholds; each
point on a curve is one iteration budget. Shaded bands are marginal 95%
bootstrap intervals; method gaps below are paired intervals.

The figures compare **pydegensac 0.3.1** with **OpenCV 5.0.0** (`RANSAC`,
`USAC_MAGSAC`) and **PoseLib 2.0.5** (with and without PROSAC) on Linux x86-64
(Intel Core i7-14700K, WSL2, gcc 13.3), single-threaded, pydegensac seed 12345.
Thresholds come from the earlier OpenCV 4.13 tuning and were not retuned for
OpenCV 5.0.

<table>
<tr><th width="50%">Fundamental matrix</th><th width="50%">Homography</th></tr>
<tr>
<td><img alt="F time-mAA, Linux x86-64" src="benchmarks/results/time_maa_f_linux.png"></td>
<td><img alt="H time-mAA, Linux x86-64" src="benchmarks/results/time_maa_h_linux.png"></td>
</tr>
</table>

**Fundamental matrix.** At each method's best budget, pydegensac reaches 0.423
mAA at 29.0 ms/pair. PoseLib PROSAC leads with 0.457 at 39.6 ms (paired gap
-0.034 mAA, 95% CI -0.055 to -0.012); OpenCV MAGSAC reaches 0.376 at 29.3 ms and
OpenCV RANSAC 0.340 at 179.0 ms. At 5,000 iterations pydegensac already reaches
0.411 at 6.3 ms.

**Homography.** On HPatchesSeq, pydegensac reaches 0.905 mAA at 1.9 ms/pair
(6,400 iterations). OpenCV MAGSAC reaches 0.927 at 0.85 ms and PoseLib PROSAC
0.930 at 7.6 ms (paired pydegensac gap -0.025 mAA, 95% CI -0.042 to -0.012),
consistent with the 2023 homography benchmark above. EVD's 8 test pairs do not
separate the methods.

**0.3.0 → 0.3.1.** Both releases built with the same toolchain, pydegensac only,
three seeds and two timing repeats per pair and budget:

| | Linux x86-64 (i7-14700K) | Apple M1 |
|---|---|---|
| F, St Peter's Square | 1.72–2.00x (1.79x summed over budgets), identical mAA | 1.18–2.10x (1.39x), identical mAA |
| H, EVD | 1.31–1.67x (1.34x), identical mAA | 1.31–1.55x (1.41x), identical mAA |
| H, HPatchesSeq | 1.20–1.24x (1.22x), mAA within 0.0002 | 1.11–1.26x (1.20x), identical mAA |

The speed-up is largest at small iteration budgets, where fixed per-call costs
dominate. On x86-64 the fundamental-matrix error and scoring loops also use AVX
when the CPU supports it, chosen at run time and bit-identical to the scalar
code (`PYDEGENSAC_SIMD=0` forces scalar); AVX alone is 1.26x on F summed over
budgets. Apple M1 runs the scalar code. On M1, pydegensac F reaches 0.425 mAA
at 38 ms/pair with 50,000 iterations (PoseLib PROSAC: 0.457 at 50 ms). On
HPatchesSeq it reaches 0.905 at 1.3 ms/pair with 6,400 iterations.
See [benchmark instructions](benchmarks/README.md) for the protocol and the
[earlier study](docs/reports/2026-08-11-ransac-benchmark.md) for historical
results and threshold-transfer caveats. Reproduce with
`cd benchmarks && python setup_data.py && ./run_current.sh linux`.

# Installation

To build and install `pydegensac`, you can use pip from Windows, macOS and Linux:

```bash
pip install pydegensac
```

> **Do not pin `pydegensac==0.1.2` if you are on numpy 2.x.** That combination
> returns **every correspondence as an inlier**, silently — on a synthetic set
> with 150 planted outliers among 300 matches it reports 300 inliers, where the
> same wheel under numpy 1.26 correctly reports 150. It is the old pybind11
> incompatibility that `0.2` was yanked for, but `0.1.2` was never yanked, so
> it is what an old pin or an unconstrained resolve on a fresh numpy can still
> land on. Use the latest release, or hold numpy below 2.0.

Or clone or download this repository and then, from within the repository, run:

```bash
pip install .
```

Building from source needs CMake, LAPACK/BLAS and a C++11 compiler (see
[Requirements](#requirements)); use the platform default compiler (Clang on
macOS). On Debian/Ubuntu the system packages are
`cmake libblas-dev liblapack-dev g++`.

To check if everything works, run the following:

```bash
cd examples
python -utt simple-example.py
```

It runs homography and fundamental matrix estimation on a sample image pair
with both OpenCV and pydegensac, and prints the inlier counts, runtimes and
estimated matrices. pydegensac runs with a fixed seed, but OpenCV's RANSAC is
unseeded and the exact numbers also depend on platform and library versions, so
no reference output is given here. pydegensac should find at least as many
inliers as OpenCV's vanilla RANSAC, in a fraction of its time.


# Example of usage

```python
import pydegensac
H, mask = pydegensac.findHomography(src_pts, dst_pts, 3.0)
F, mask = pydegensac.findFundamentalMatrix(src_pts, dst_pts, 3.0)

```

See also this [notebook](examples/simple-example.ipynb) with simple example

And this [notebook](examples/how-to-use-detailed.ipynb) with detailed explanation of possible options


# Requirements

- Python 3.8 or higher
- CMake 3.9 or higher
- LAPACK
- BLAS (OpenBLAS, MKL, Atlas, ...)
- A modern compiler with C++11 support


## Citation

Please cite us if you use this code:

    @InProceedings{Chum2003,
    author="Chum, Ond{\v{r}}ej and Matas, Ji{\v{r}}{\'i} and Kittler, Josef",
    title="Locally Optimized RANSAC",
    booktitle="Pattern Recognition",
    year="2003",
    }
    
    @inproceedings{Chum2005,
    author = {Chum, Ondrej and Werner, Tomas and Matas, Jiri},
    title = {Two-View Geometry Estimation Unaffected by a Dominant Plane},
    booktitle = {CVPR},
    year = {2005},
    }
    
    @article{Mishkin2015MODS,
          title = "MODS: Fast and robust method for two-view matching ",
          journal = "Computer Vision and Image Understanding ",
          year = "2015",
          issn = "1077-3142",
          doi = "http://dx.doi.org/10.1016/j.cviu.2015.08.005",
          url = "http://www.sciencedirect.com/science/article/pii/S1077314215001800",
          author = "Dmytro Mishkin and Jiri Matas and Michal Perdoch"
    }
    


    
# Acknowledgements

This wrapper part is based on great [Benjamin Jack `python_cpp_example`](https://github.com/benjaminjack/python_cpp_example).
