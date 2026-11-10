"""Compile and run the guarded 3x3 rank-two projection regression harness."""
from pathlib import Path
import subprocess

import pytest


def test_rank2_projection_numerics_and_fallbacks(tmp_path):
    ccompiler = pytest.importorskip("setuptools._distutils.ccompiler")
    sysconfig = pytest.importorskip("setuptools._distutils.sysconfig")
    compiler = ccompiler.new_compiler()
    sysconfig.customize_compiler(compiler)
    root = Path(__file__).resolve().parents[1]
    source = root / "src/pydegensac/degensac"
    objects = compiler.compile(
        [str(root / "tests/c/test_rank2_projection.c"), str(source / "rank2.c")],
        output_dir=str(tmp_path), include_dirs=[str(source)],
    )
    executable = str(tmp_path / "test_rank2_projection")
    libraries = [] if compiler.compiler_type == "msvc" else ["m"]
    compiler.link_executable(objects, executable, libraries=libraries)
    subprocess.run([compiler.executable_filename(executable)], check=True)
