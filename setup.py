
from setuptools import setup, Extension
from Cython.Build import cythonize
import os
import numpy

#########################
### INSTALLATION NOTE ###
#########################

# PHAT is often *very difficult* to install. Before installing, please ensure
# that the following are installed and accessible in your system:
#
#   • PyBind11, setuptools, wheel (via pip)
#   • g++ (preferably >12)
#
# When installing PHAT, please use the arguments
#
#   pip install --use-deprecated=legacy-resolver --no-binary :all: phat
#
# To ensure installation.

if not "CC" in os.environ: os.environ["CC"] = "gcc-14"
if not "CXX" in os.environ: os.environ["CXX"] = "g++-14"

extensions = [
	Extension(
		"*",
		["ateams/**/*.pyx"],
		include_dirs=[numpy.get_include()],
		extra_compile_args=["-O4"], # just C
		# extra_compile_args=["-std=c++20", "-O3"] # C++
	)
]

setup(
    ext_modules=cythonize(extensions, annotate=True)
)