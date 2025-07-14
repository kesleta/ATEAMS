
from ..common cimport Table, Lookup, BoundaryMatrix, Index, Vector, Set

cdef extern from "LinBoxMethods.h":
	Vector LanczosKernelSample(Index coboundary, int M, int N, int p, int maxTries) except +
	Set ComputePercolationEvents(Table addition, Table multiplication, Lookup negation, Lookup inversion, BoundaryMatrix Boundary, Index breaks, int cellCount, int topDimension, int homology) noexcept
