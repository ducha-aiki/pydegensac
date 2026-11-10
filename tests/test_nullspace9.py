"""Compile and run the exact fixed-size 9x9 nullspace regression harness."""
from pathlib import Path
import subprocess

import pytest


def test_nullspace9_matches_generic_solver(tmp_path):
    ccompiler = pytest.importorskip("setuptools._distutils.ccompiler")
    sysconfig = pytest.importorskip("setuptools._distutils.sysconfig")
    compiler = ccompiler.new_compiler()
    sysconfig.customize_compiler(compiler)
    root = Path(__file__).resolve().parents[1]
    source = root / "src/pydegensac/degensac"
    objects = compiler.compile(
        [str(root / "tests/c/test_nullspace9.c"), str(source / "utools.c")],
        output_dir=str(tmp_path), include_dirs=[str(source)], macros=[("NDEBUG",)],
    )
    executable = str(tmp_path / "test_nullspace9")
    libraries = [] if compiler.compiler_type == "msvc" else ["m"]
    compiler.link_executable(objects, executable, libraries=libraries)
    subprocess.run([compiler.executable_filename(executable)], check=True)
