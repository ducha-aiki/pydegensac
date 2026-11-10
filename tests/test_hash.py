"""Exercise hash collisions, overflow storage, and repeated cleanup."""
from pathlib import Path
import subprocess

import pytest


def test_hash_storage_and_lookup(tmp_path):
    ccompiler = pytest.importorskip("setuptools._distutils.ccompiler")
    sysconfig = pytest.importorskip("setuptools._distutils.sysconfig")
    compiler = ccompiler.new_compiler()
    sysconfig.customize_compiler(compiler)
    root = Path(__file__).resolve().parents[1]
    source_dir = root / "src/pydegensac/degensac"
    objects = compiler.compile(
        [str(root / "tests/c/test_hash.c"), str(source_dir / "hash.c")],
        output_dir=str(tmp_path), include_dirs=[str(source_dir)],
        # A one-element tuple undefines NDEBUG, keeping C assertions active.
        macros=[("NDEBUG",)],
    )
    executable = str(tmp_path / "test_hash")
    compiler.link_executable(objects, executable)
    subprocess.run([compiler.executable_filename(executable)], check=True)
