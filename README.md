# pydegensac

This repository contains an Python wrapper of RANSAC for homography and fundamental matrix estimation
from sparse correspondences. It implements [LO-RANSAC](https://link.springer.com/chapter/10.1007/978-3-540-45243-0_31) and [DEGENSAC](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.466.2719&rep=rep1&type=pdf).

It was originally located in [https://github.com/ducha-aiki/pyransac](https://github.com/ducha-aiki/pyransac), but was renamed to avoid conflict with already existing [pyransac](https://pypi.org/project/pyransac/) in pypi from other author.

# Performance

Vanilla pydegensac implementation is marginally better than OpenCV one and with degeneracy-check enabled (DEGENSAC) it is the state of the art,
according to the recent study Yin et.al."[Image Matching across Wide Baselines: From Paper to Practice](https://arxiv.org/abs/2003.01587.pdf)", 2020.

![IMW-benchmark](img/ransacs.png)


![IMW-Challenge](img/ransacs2.png)


# Installation

To build and install `pydegensac`, you can use pip from Windows, macOS and Linux:

```bash
pip install pydegensac
```

Or clone or download this repository and then, from within the repository, run:

```bash
python3 ./setup.py install
```

or

```bash
pip3 install .
```

To check if everything works, run the following:

```bash
cd examples
python -utt simple-example.py
```

You should see the following output:

```
Running homography estimation
cv2 found 40 inliers
OpenCV runtime 0.02355  sec
pydegensac found 78 inliers
pydegensac runtime 0.00320  sec
H =  [[ 5.59934334e-03 -2.36037104e-03 -2.78369679e+01]
 [ 4.86321171e-02 -1.24542142e-01 -1.00600649e+01]
 [ 1.95536148e-04  9.43300063e-06 -1.76685691e-01]]
Running fundamental matrix estimation
cv2 found 32 inliers
OpenCV runtime 0.67554  sec
pydegensac found 44 inliers
pydegensac 0.04702  sec
F =  [[-7.35044984e-04 -2.72572333e-03  1.38155992e+00]
 [ 1.43946998e-03  2.33120834e-05 -7.88961637e-01]
 [-3.35556093e-01  1.00000000e+00 -1.78675406e+02]]
```

# Building hints from Tomasz Malisiewicz

1. Compiling pydegensac without a system-wide install.

```bash
python3 ./setup.py build
```

2. Compiling on Mac OS X computer
Use GCC instead of Clang. The most recent version on my machine (installed via brew) is gcc-8. Try this:

```bash
CC=gcc-8 python3 ./setup.py build
```

3. Compiling on Ubuntu 18.04
You need LAPACK and a few other libraries and I always forget those specific package names. Take a look at my pydegensac Dockerfile to see the exact packages you need to apt install on an Ubuntu 18.04 system (https://github.com/quantombone/pydegensac-dockerfile/blob/master/Dockerfile)

```bash
FROM ubuntu:18.04
```

## update system
```bash
RUN apt-get clean
RUN apt-get update
RUN apt-get install -qy \
    git python3 python3-setuptools python3-dev
RUN apt-get install -y cmake libblas-dev liblapack-dev gfortran
RUN apt-get install -y g++ gcc
```

## download and build pydegensac
```
RUN git clone https://github.com/ducha-aiki/pydegensac.git
WORKDIR pydegensac
RUN python3 ./setup.py build
```

## copy built assets into target directory (which will be a -v volume)
```docker
CMD cp -R /pydegensac/build/lib.linux-x86_64-3.6/pydegensac /target_directory
```

# dockerfile

https://github.com/quantombone/pydegensac-dockerfile


# Example of usage

```python
import pydegensac
H, mask = pydegensac.findHomography(src_pts, dst_pts, 3.0)
F, mask = pydegensac.findFundamentalMatrix(src_pts, dst_pts, 3.0)

```

See also this [notebook](examples/simple-example.ipynb) with simple example

And this [notebook](examples/how-to-use-detailed.ipynb) with detailed explanation of possible options


# Requirements

- Python 3
- CMake 2.8.12 or higher
- LAPACK, 
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
