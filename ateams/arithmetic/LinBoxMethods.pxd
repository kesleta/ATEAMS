
from ..common cimport Table, Lookup, BoundaryMatrix, Index, Set

cdef extern from "LinBoxMethods.h":
	Index ReducedKernelSample(Index coboundary, int M, int N, int p) except +
	Set ComputePercolationEvents(Table addition, Table multiplication, Lookup negation, Lookup inversion, BoundaryMatrix Boundary, Index breaks, int cellCount) noexcept
	Set LinearComputePercolationEvents(int field, Lookup addition, Lookup multiplication, Lookup negation, Lookup inversion, BoundaryMatrix Boundary, Index breaks, int cellCount, int dimension) noexcept
	Set ZpComputePercolationEvents(int field, BoundaryMatrix Boundary, Index breaks, int cellCount) noexcept
