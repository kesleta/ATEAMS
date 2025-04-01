
# cython: language_level=3str, initializedcheck=False, c_api_binop_methods=True, nonecheck=False, profile=True, cdivision=True, boundscheck=False, wraparound=False
# define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
# distutils: define_macros=CYTHON_TRACE_noexcept nogil=1
# distutils: language=c++

from .common cimport FFINT, FLAT, TABLE

from libcpp cimport bool
from libcpp.vector cimport vector as Vector
from libcpp.unordered_set cimport unordered_set as Set
from libcpp.unordered_map cimport unordered_map as Map


cdef class Matrix:
	cdef TABLE addition
	cdef FLAT negation
	cdef TABLE multiplication
	cdef FLAT inverses
	cdef TABLE data
	cdef bool parallel
	cdef str schedule
	cdef int cores
	cdef int minBlockSize
	cdef int maxBlockSize
	cdef int blockSize

	cdef Map[int, Set[int]] columns
	cdef Vector[int] shape
	cdef Vector[Vector[int]] blockSchema

	cdef void _initializeColumns(self) 
	cdef TABLE ToArray(self) 

	cdef void SwapRows(self, int i, int j) noexcept nogil
	cdef void AddRows(self, int i, int j, int MINCOL, int MAXCOL, FFINT ratio) noexcept nogil
	cdef void BlockAddRows(self, int i, int j, FFINT ratio, int start, int stop, str schedule) noexcept nogil
	cdef void MultiplyRow(self, int i, FFINT q) noexcept nogil
	cdef void RescanColumns(self, int MINCOL) noexcept nogil
	cdef int PivotRow(self, int c, int pivots) noexcept nogil
	cdef int HighestZeroRow(self, int AUGMENT=*) 

	cdef void RREF(self, int AUGMENT=*) noexcept nogil
