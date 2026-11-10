"""Compile and run the exact seven-point padded-nullspace regression harness."""
from pathlib import Path
import subprocess

import pytest


def test_nullspace7x9_matches_generic_padded_solver(tmp_path):
    ccompiler = pytest.importorskip("setuptools._distutils.ccompiler")
    sysconfig = pytest.importorskip("setuptools._distutils.sysconfig")
    compiler = ccompiler.new_compiler()
    sysconfig.customize_compiler(compiler)
    root = Path(__file__).resolve().parents[1]
    source = root / "src/pydegensac/degensac"
    objects = compiler.compile(
        [str(root / "tests/c/test_nullspace7x9.c"), str(source / "utools.c")],
        output_dir=str(tmp_path), include_dirs=[str(source)],
        macros=[("NDEBUG",)],
    )
    executable = str(tmp_path / "test_nullspace7x9")
    libraries = [] if compiler.compiler_type == "msvc" else ["m"]
    compiler.link_executable(objects, executable, libraries=libraries)
    subprocess.run([compiler.executable_filename(executable)], check=True)
