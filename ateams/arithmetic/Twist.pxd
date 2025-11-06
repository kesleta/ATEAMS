
# distutils: language=c++

from ..common cimport INDEXFLAT, BoundaryMatrix, Index, Set, Table, Lookup, PersistencePairs



# distutils: language=c++

from ..common cimport INDEXFLAT, BoundaryMatrix, Index, Set, Table, Lookup, Bases, Basis, bool


cdef class Twist:
	cdef int characteristic, cellCount, dimension, topDimension
	cdef bool __DEBUG
	cdef Table addition, multiplication
	cdef Lookup negation, inversion, flatAddition, flatMultiplication

	cdef Index fullBoundary, breaks
	cdef BoundaryMatrix referenceBoundary, workingBoundary, partialBoundary, partialCoboundary, augmentedCoboundary
	cdef Bases bases
	cdef Basis cobasis

	cdef void __arithmetic(self) noexcept
	cdef BoundaryMatrix __transpose(self, BoundaryMatrix A, int columns) noexcept
	cdef BoundaryMatrix __asRows(self, BoundaryMatrix A, int rows) noexcept
	cdef BoundaryMatrix __asColumns(self, BoundaryMatrix A, int columns) noexcept

	cdef BoundaryMatrix FillBoundaryMatrix(self, INDEXFLAT boundary) noexcept
	cdef BoundaryMatrix ReindexBoundaryMatrix(self, INDEXFLAT filtration) noexcept
	cdef BoundaryMatrix PartialBoundaryMatrix(self, int dimension) noexcept
	cdef Index ReindexPartialFiltration(self, INDEXFLAT filtration) noexcept
	cdef BoundaryMatrix ReindexPartialBoundaryMatrix(self, BoundaryMatrix &boundary, INDEXFLAT partial) noexcept
	
	cpdef Set LinearComputePercolationEvents(self, INDEXFLAT filtration) noexcept
	cpdef Set RankComputePercolationEvents(self, INDEXFLAT filtration) noexcept
	cpdef Set SolveComputePercolationEvents(self, INDEXFLAT filtration) noexcept

	cpdef Basis LinearComputeBasis(self) noexcept
	cpdef Basis LinearComputeCobasis(self) noexcept



cpdef PersistencePairs PHATComputePersistencePairs(INDEXFLAT boundary, INDEXFLAT filtration, int homology, INDEXFLAT breaks) noexcept
