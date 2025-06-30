
from setuptools import setup, Extension
from Cython.Build import cythonize
from Cython.Compiler.Options import get_directive_defaults as GDD
from Cython.Compiler import Options
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

# try:
# 	if os.environ["CC"]: os.environ["CC"] = "gcc-14"
# 	if os.environ["CXX"]: os.environ["CXX"] = "g++-14"
# except:
# 	os.environ["CC"] = "gcc-12"
# 	os.environ["CXX"] = "g++-12"

# os.environ["CC"] = "gcc"
# os.environ["CXX"] = "g++"

DD = GDD()
DD["linetrace"] = True
DD["binding"] = True

extensions = [
	Extension(
		"*",
		["ateams/**/*.pyx"],
		include_dirs=[numpy.get_include()],
		extra_compile_args=["-std=c++20", "-O3", "-lz"],
		define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
		language="c++",
		libraries=["FastSample"]
	)
]

## STORE THE SHARED LIBRARY IN /usr/local/lib FOR EASE-OF-USE

setup(
    ext_modules=cythonize(
		extensions,
		annotate=True,
		language_level="3",
		compiler_directives=dict(
			initializedcheck=False,
			c_api_binop_methods=True,
			nonecheck=False,
			profile=True,
			cdivision=True,
			binding=True,
			linetrace=True,
			boundscheck=False,
			wraparound=False
		)
	)
)