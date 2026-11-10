"""Compile and run the exact eight-point nullspace regression harness."""
from pathlib import Path
import subprocess

import pytest


def test_svdnull8_matches_full_svd(tmp_path):
    ccompiler = pytest.importorskip("setuptools._distutils.ccompiler")
    sysconfig = pytest.importorskip("setuptools._distutils.sysconfig")
    compiler = ccompiler.new_compiler()
    sysconfig.customize_compiler(compiler)
    root = Path(__file__).resolve().parents[1]
    matutls = root / "src/pydegensac/matutls"
    sources = [
        root / "tests/c/test_svdnull8.c",
        matutls / "svduv.c",
        matutls / "ldumat.c",
        matutls / "ldvmat.c",
        matutls / "qrbdv.c",
    ]
    objects = compiler.compile(
        [str(source) for source in sources],
        output_dir=str(tmp_path),
        include_dirs=[str(matutls)],
        # A one-element tuple undefines NDEBUG, keeping C assertions active.
        macros=[("NDEBUG",)],
    )
    executable = str(tmp_path / "test_svdnull8")
    # MSVC has no separate libm; its math functions live in the CRT.
    libraries = [] if compiler.compiler_type == "msvc" else ["m"]
    compiler.link_executable(objects, executable, libraries=libraries)
    subprocess.run([compiler.executable_filename(executable)], check=True)
