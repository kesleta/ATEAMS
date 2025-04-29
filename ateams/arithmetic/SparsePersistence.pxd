
from .common cimport FFINT, FLATCONTIG, TABLECONTIG, INDEXTABLE, INDEXFLAT

import numpy as np
cimport numpy as np

from libcpp cimport bool
from libcpp.vector cimport vector as Vector
from libcpp.unordered_set cimport unordered_set as Set
from libcpp.set cimport set as OrderedSet
from libcpp.unordered_map cimport unordered_map as Map


cdef class Persistence:
	cdef TABLECONTIG addition
	cdef TABLECONTIG subtraction
	cdef TABLECONTIG multiplication
	cdef FLATCONTIG inverse
	cdef FLATCONTIG negation
	cdef FFINT characteristic

	cdef void __arithmetic(self) noexcept
	cdef void __flushDataStructures(self) noexcept

	cdef INDEXTABLE tranches
	cdef INDEXFLAT dimensions
	cdef int homology
	cdef int cellCount
	cdef int vertexCount
	cdef int higherCellCount
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

	cdef Vector[Vector[int]] boundary
	cdef Vector[Vector[int]] _boundary
	cdef Vector[int] _dimensions
	cdef Vector[Vector[int]] _tranches

	cdef Vector[int] markedIterable
	cdef Set[int] marked
	cdef Vector[int] premarked
	
	cpdef Vector[Vector[int]] ReindexBoundary(self, INDEXFLAT filtration) noexcept
	cdef Vector[Vector[int]] ReindexSubBoundary(self, INDEXFLAT subcomplex) noexcept
	cdef Vector[Vector[int]] Vectorize(self, list[list[int]] flattened) noexcept

	cdef OrderedSet[int] RemoveUnmarkedCells(self, int cell, OrderedSet[int] faces, Map[int,FFINT] &faceCoefficients) noexcept
	cdef OrderedSet[int] Eliminate(self, int youngest, OrderedSet[int] faces, Map[int,FFINT] &faceCoefficients) noexcept
	cdef OrderedSet[int] ReducePivotRow(self, int cell, OrderedSet[int] faces, Map[int,FFINT] &faceCoefficients) noexcept
	cpdef OrderedSet[int] ComputePercolationEvents(self, INDEXFLAT filtration) noexcept
	cpdef Vector[int] ComputeBettiNumbers(self, INDEXFLAT subcomplex) noexcept
