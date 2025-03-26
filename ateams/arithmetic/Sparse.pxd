
# cython: language_level=3str, initializedcheck=False, c_api_binop_methods=True, nonecheck=False, profile=True, cdivision=True, wraparound=False, boundscheck=False
# define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
# distutils: define_macros=CYTHON_TRACE_NOGIL=1
# distutils: language=c++

from .common cimport FFINT, FLAT, TABLE

from libcpp.vector cimport vector as Vector
from libcpp.set cimport set as Set
from libcpp.map cimport map as Map


cdef class Matrix:
	cdef TABLE addition
	cdef FLAT negation
	cdef TABLE multiplication
	cdef FLAT inverses
	cdef TABLE _original

	cdef Map[int, Map[int, FFINT]] rows
	cdef Map[int, Set[int]] columns
	cdef Vector[int] shape

	cdef void _initializeRows(self, TABLE A) noexcept
	cdef TABLE ToArray(self) noexcept nogil

	cdef void SwapRows(self, int i, int j) noexcept nogil
	cdef void AddRows(self, int i, int j, FFINT ratio) noexcept nogil
	cdef void MultiplyRow(self, int i, FFINT q) noexcept nogil
	cdef int PivotRow(self, int c, int pivots) noexcept nogil
	cdef int HighestZeroRow(self, int AUGMENT=*) noexcept

	cdef void RREF(self, int AUGMENT=*) noexcept
