import os
import re
import sys
import platform
import subprocess
import multiprocessing
import shutil
import glob

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: " +
                ", ".join(e.name for e in self.extensions))

        if self.is_windows():
            cmake_version = LooseVersion(
                re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    @staticmethod
    def is_windows():
        tag = platform.system().lower()
        return tag == "windows"

    @staticmethod
    def is_linux():
        tag = platform.system().lower()
        return tag == "linux"

    def build_extension(self, ext):
        extdir = os.path.abspath(
            os.path.dirname(self.get_ext_fullpath(ext.name)))
        if not os.path.isdir(extdir):
            os.makedirs(extdir)

        cfg = 'Release'
        build_args = ['--config', cfg]

        cpu_count = max(2, multiprocessing.cpu_count() // 2)
        build_args += ['--', '-j{0}'.format(cpu_count)]

        # contrib folder
        contrib_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        # build smt-switch
        contrib_smt_switch = [os.path.join(contrib_path, "setup-smt-switch.sh"), "--python"]
        subprocess.check_call(contrib_smt_switch)
        # build btor2
        contrib_btor2 = os.path.join(contrib_path, "setup-btor2tools.sh")
        subprocess.check_call(contrib_btor2)
        # build bison
        contrib_bison = os.path.join(contrib_path, "setup-bison.sh")
        subprocess.check_call(contrib_bison)
        # build flex
        contrib_flex = os.path.join(contrib_path, "setup-flex.sh")
        subprocess.check_call(contrib_flex)
        # build coreir
        contrib_coreir = os.path.join(contrib_path, "setup-coreir.sh")
        subprocess.check_call(contrib_coreir)

        # configure
        root_dir = os.path.dirname(contrib_path)
        build_dir = os.path.join(root_dir, "build")
        # to avoid multiple build, only call reconfigure if we couldn't find the makefile
        # for python
        python_make_dir = os.path.join(build_dir, "python")
        if not os.path.isfile(os.path.join(python_make_dir, "Makefile")):
            configure_path = os.path.join(root_dir, "configure.sh")
            configure_args = [configure_path, "--python", "--with-coreir"]
            subprocess.check_call(configure_args, cwd=root_dir)

        # build the main library
        subprocess.check_call(
            ['cmake', '--build', '.', "--target", "pono"] + build_args, cwd=build_dir)
        # build the python binding
        python_build_dir = os.path.join(build_dir, "python")
        subprocess.check_call(["make"], cwd=python_build_dir)
        # the build folder gets cleaned during the config, need to create it again
        # this is necessary since "build" is a python dist folder
        if not os.path.isdir(extdir):
            os.mkdir(extdir)

        for lib_filename in glob.glob(os.path.join(python_build_dir, "pono.*")):
            if os.path.splitext(lib_filename)[1] == ".cxx":
                continue
            dst_filename = os.path.join(extdir, os.path.basename(lib_filename))
            shutil.copy(lib_filename, dst_filename)


setup(
    name='pono',
    version='0.1.1',
    author='Makai Mann',
    ext_modules=[CMakeExtension('pono')],
    cmdclass=dict(build_ext=CMakeBuild),
    long_description='python bindings for Pono (next generation of CoSA)',
    url='https://github.com/upscale-project/pono',
    license='BSD',
    install_requires=['smt-switch'],
    tests_require=['pytest'],
    zip_safe=False,
)

