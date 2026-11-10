"""Compile and run the exact index-only scorer regression harness."""
from pathlib import Path
import subprocess

import pytest


def test_inlidxs_only_matches_full_scorer(tmp_path):
    ccompiler = pytest.importorskip("setuptools._distutils.ccompiler")
    sysconfig = pytest.importorskip("setuptools._distutils.sysconfig")
    compiler = ccompiler.new_compiler()
    sysconfig.customize_compiler(compiler)
    root = Path(__file__).resolve().parents[1]
    source = root / "src/pydegensac/degensac"
    objects = compiler.compile(
        [str(root / "tests/c/test_inlidxs_only.c"), str(source / "rtools.c"),
         str(source / "bsd_random.c")],
        output_dir=str(tmp_path), include_dirs=[str(source)],
        macros=[("NDEBUG",)],
    )
    executable = str(tmp_path / "test_inlidxs_only")
    libraries = [] if compiler.compiler_type == "msvc" else ["m"]
    compiler.link_executable(objects, executable, libraries=libraries)
    subprocess.run([compiler.executable_filename(executable)], check=True)
