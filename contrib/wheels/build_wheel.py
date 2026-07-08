#!/usr/bin/env python3
import multiprocessing
import platform
import re
import shutil
import subprocess
from distutils.version import LooseVersion
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name: str) -> None:
        Extension.__init__(self, name, sources=[])


class CMakeBuild(build_ext):
    def run(self) -> None:
        cmake_path = shutil.which("cmake")
        if cmake_path is None:
            msg = "cmake not found"
            raise RuntimeError(msg)
        self.cmake_path = cmake_path
        try:
            out = subprocess.check_output([self.cmake_path, "--version"])
        except OSError as err:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            ) from err

        if self.is_windows():
            match = re.search(r"version\s*([\d.]+)", out.decode())
            if match is None:
                msg = f"Unexpected CMake version string: {out.decode()}"
                raise RuntimeError(msg)
            cmake_version = LooseVersion(match.group(1))
            if cmake_version < "3.1.0":
                msg = "CMake >= 3.1.0 is required on Windows"
                raise RuntimeError(msg)

        for ext in self.extensions:
            self.build_extension(ext)

    @staticmethod
    def is_windows() -> bool:
        tag = platform.system().lower()
        return tag == "windows"

    @staticmethod
    def is_linux() -> bool:
        tag = platform.system().lower()
        return tag == "linux"

    def build_extension(self, ext: Extension) -> None:
        extdir = Path(self.get_ext_fullpath(ext.name)).parent.resolve()
        if not extdir.is_dir():
            extdir.mkdir(parents=True)

        cfg = "Release"
        build_args = ["--config", cfg]

        cpu_count = max(2, multiprocessing.cpu_count() // 2)
        build_args += ["--", f"-j{cpu_count}"]

        # contrib folder
        contrib_path = Path(__file__).resolve().parent.parent
        # build smt-switch
        subprocess.check_call(
            [
                contrib_path / "setup-smt-switch.sh",
                "--python",
            ]
        )
        # build btor2
        subprocess.check_call(contrib_path / "setup-btor2tools.sh")
        # build bison
        subprocess.check_call(contrib_path / "setup-bison.sh")
        # build flex
        subprocess.check_call(contrib_path / "setup-flex.sh")
        # build coreir
        subprocess.check_call(contrib_path / "setup-coreir.sh")

        # configure
        root_dir = contrib_path.parent
        build_dir = root_dir / "build"
        # to avoid multiple builds, only reconfigure if we couldn't find the makefile
        # for python
        python_build_dir = build_dir / "python"
        python_makefile = python_build_dir / "Makefile"
        if not python_makefile.is_file():
            configure_path = root_dir / "configure.sh"
            subprocess.check_call(
                [configure_path, "--python", "--with-coreir"], cwd=root_dir
            )

        # build the main library
        subprocess.check_call(
            [self.cmake_path, "--build", ".", "--target", "pono", *build_args],
            cwd=build_dir,
        )
        # build the python binding
        make_path = shutil.which("make")
        if make_path is None:
            msg = "make not found"
            raise RuntimeError(msg)
        subprocess.check_call([make_path], cwd=python_build_dir)
        # the build folder gets cleaned during the config, need to create it again
        # this is necessary since "build" is a python dist folder
        if not extdir.is_dir():
            extdir.mkdir()

        for lib_filename in python_build_dir.glob("pono.*"):
            if lib_filename.suffix == ".cxx":
                continue
            dst_filename = extdir / lib_filename.name
            shutil.copy(lib_filename, dst_filename)


setup(
    name="pono",
    version="0.1.1",
    author="Makai Mann",
    ext_modules=[CMakeExtension("pono")],
    cmdclass={"build_ext": CMakeBuild},
    long_description="python bindings for Pono (next generation of CoSA)",
    url="https://github.com/upscale-project/pono",
    license="BSD",
    install_requires=["smt-switch"],
    tests_require=["pytest"],
    zip_safe=False,
)
