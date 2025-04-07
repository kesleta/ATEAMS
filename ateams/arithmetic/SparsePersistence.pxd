
# cython: language_level=3str, initializedcheck=False, c_api_binop_methods=True, nonecheck=False, profile=True, cdivision=True
# cython: boundscheck=False, wraparound=False
# cython: binding=True, linetrace=True
# define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
# distutils: language=c++

from .common cimport FFINT, FLAT, TABLE, INDEXTABLE, INDEXFLAT

import numpy as np
cimport numpy as np

from libcpp cimport bool
from libcpp.vector cimport vector as Vector
from libcpp.unordered_set cimport unordered_set as Set
from libcpp.set cimport set as OrderedSet
from libcpp.unordered_map cimport unordered_map as Map


cdef class Persistence:
	cdef TABLE addition
	cdef TABLE subtraction
	cdef TABLE multiplication
	cdef FLAT inverse
	cdef FLAT negation
	cdef FFINT characteristic

	cdef void __arithmetic(self) noexcept
	cdef void __flushDataStructures(self) noexcept

	cdef INDEXTABLE tranches
	cdef INDEXFLAT dimensions
	cdef int homology
	cdef int cellCount
	cdef int vertexCount
	cdef int defaultRowSize
	cdef int tagged
	cdef int low
	cdef int high

	cdef int cores
	cdef int minBlockSize
	cdef int maxBlockSize
	cdef bool parallel

	cdef Vector[OrderedSet[int]] columnEntries
	cdef Vector[Vector[int]] columnEntriesIterable
	cdef Vector[Map[int,FFINT]] columnEntriesCoefficients

	cdef Vector[Vector[int]] boundary;
	cdef Vector[int] markedIterable;
	cdef Set[int] marked;
	cdef Vector[int] premarked;

	cdef Vector[Vector[int]] Vectorize(self, list[list[int]] flattened) noexcept

	cdef OrderedSet[int] RemoveUnmarkedCells(self, int cell, OrderedSet[int] faces, Map[int,FFINT] &faceCoefficients) noexcept
	cdef OrderedSet[int] Eliminate(self, int youngest, OrderedSet[int] faces, Map[int,FFINT] &faceCoefficients) noexcept
	cdef OrderedSet[int] ReducePivotRow(self, int cell, OrderedSet[int] faces, Map[int,FFINT] &faceCoefficients) noexcept
	cpdef OrderedSet[int] ComputePercolationEvents(self, INDEXFLAT filtration, list[list[int]] flattened) noexcept

