
# cython: language_level=3str, initializedcheck=False, c_api_binop_methods=True, nonecheck=False, profile=True, cdivision=True
# cython: boundscheck=False, wraparound=False
# cython: binding=True, linetrace=True
# define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
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
	cdef int zero
	cdef int cores
	cdef int minBlockSize
	cdef int maxBlockSize
	cdef int blockSize

	cdef Map[int, Set[int]] columns
	cdef Map[int, Map[int, Set[int]]] positiveThreadCache
	cdef Map[int, Map[int, Set[int]]] negativeThreadCache
	cdef Vector[int] shape
	cdef Vector[Vector[int]] blockSchema

	cdef void _initializeColumns(self) noexcept
	cdef void _initializeThreadCaches(self) noexcept
	cdef void _flushThreadCaches(self, int MINROW, int MAXROW) noexcept
	cdef TABLE ToArray(self) noexcept

	cdef void SwapRows(self, int i, int j) noexcept
	cdef void AddRows(self, int i, int j, int MINCOL, int MAXCOL, FFINT ratio) noexcept
	cdef void ParallelAddRows(self, int i, int j, int start, int stop, FFINT ratio) noexcept nogil
	cdef void MultiplyRow(self, int i, FFINT q) noexcept
	cdef Vector[Vector[int]] recomputeBlockSchema(self, int start, int stop) noexcept
	cdef int PivotRow(self, int c, int pivots) noexcept
	cdef int HighestZeroRow(self, int AUGMENT=*) noexcept

	cdef void RREF(self, int AUGMENT=*) noexcept
