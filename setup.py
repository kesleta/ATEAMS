
from setuptools import setup, Extension
from Cython.Build import cythonize
from Cython.Compiler.Options import get_directive_defaults as GDD
import numpy

DD = GDD()
DD["linetrace"] = True
DD["binding"] = True

## STORE THE SHARED LIBRARY IN /usr/local/lib FOR EASE-OF-USE
LIBRARY_PATH = "/usr/local/lib"
INCLUDE_PATH = "/usr/local/include"

extensions = [
	Extension(
		"*",
		["ateams/**/*.pyx"],
		include_dirs=[numpy.get_include(), INCLUDE_PATH],
		library_dirs=[LIBRARY_PATH],
		extra_compile_args=["-std=c++20", "-O2"],
		define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
		language="c++",
		libraries=["LinBoxMethods", "PHATMethods", "spasm"]
	)
]

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