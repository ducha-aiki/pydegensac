#! /usr/bin/env python3
import os
import re
import sys
import platform
import subprocess

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        super().__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError as e:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(ext.name for ext in self.extensions)
            ) from e

        if platform.system() == "Windows":
            m = re.search(r"version\s*([\d.]+)", out.decode())
            ver = tuple(int(x) for x in m.group(1).split(".")[:3]) if m else (0,)
            if ver < (3, 1):
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))

        cfg = "Debug" if self.debug else "Release"
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",
        ]

        # --- macOS arch handling (cibuildwheel-friendly) ---
        # cibuildwheel sets _PYTHON_HOST_PLATFORM per build config, e.g.:
        #   macosx-11.0-x86_64
        #   macosx-11.0-arm64
        # If CMAKE_OSX_ARCHITECTURES is not set, derive it from that.
        if platform.system() == "Darwin" and "CMAKE_OSX_ARCHITECTURES" not in os.environ:
            host = os.environ.get("_PYTHON_HOST_PLATFORM")
            if host:
                # last component is usually the arch
                os.environ["CMAKE_OSX_ARCHITECTURES"] = host.split("-")[-1]

        osx_archs = os.environ.get("CMAKE_OSX_ARCHITECTURES")
        if platform.system() == "Darwin" and osx_archs:
            # Pass explicitly to CMake (environment-only is not reliably honored)
            cmake_args.append(f"-DCMAKE_OSX_ARCHITECTURES={osx_archs}")

        build_args = ["--config", cfg]

        if platform.system() == "Windows":
            # Multi-config generators (VS)
            cmake_args += [f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{cfg.upper()}={extdir}"]
            if sys.maxsize > 2**32:
                cmake_args += ["-A", "x64"]
            build_args += ["--", "/m"]
        else:
            # Make / Ninja
            build_args += ["--", "-j2"]

        env = os.environ.copy()
        env["CXXFLAGS"] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get("CXXFLAGS", ""), self.distribution.get_version()
        )

        # Build directory: isolate per-arch on macOS to prevent cache contamination
        build_temp = self.build_temp
        if platform.system() == "Darwin" and osx_archs:
            safe_arch = osx_archs.replace(";", "_").replace(" ", "_")
            build_temp = f"{build_temp}-{safe_arch}"

        if not os.path.exists(build_temp):
            os.makedirs(build_temp)

        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=build_temp, env=env)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=build_temp)

        print()  # Add an empty line for cleaner output


requirements = ["numpy"]

from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="pydegensac",
    version="0.3.0",
    author="Ondra Chum, Dmytro Mishkin",
    author_email="ducha.aiki@gmail.com",
    license="MIT",
    url="https://github.com/ducha-aiki/pydegensac",
    description="Advanced RANSAC (DEGENSAC) with bells and whistles for H and F estimation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages("src"),
    package_dir={"": "src"},
    ext_modules=[CMakeExtension("pydegensac/pydegensac")],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
    install_requires=requirements,
)
