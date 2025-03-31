
from setuptools import setup, Extension
from Cython.Build import cythonize
from Cython.Compiler.Options import get_directive_defaults as GDD
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
#

try:
	if os.environ["CC"]: os.environ["CC"] = "gcc-14"
	if os.environ["CXX"]: os.environ["CXX"] = "g++-14"
except:
	os.environ["CC"] = "gcc-12"
	os.environ["CXX"] = "g++-12"

DD = GDD()
DD["linetrace"] = True
DD["binding"] = True

extensions = [
	Extension(
		"*",
		["ateams/**/*.pyx"],
		include_dirs=[numpy.get_include(), "ateams/arithmetic"],
		extra_compile_args=["-std=c++20", "-O4", "-fopenmp"],
        extra_link_args=["-fopenmp"],
	)
]

setup(
    ext_modules=cythonize(extensions, annotate=True, language_level="3")
)