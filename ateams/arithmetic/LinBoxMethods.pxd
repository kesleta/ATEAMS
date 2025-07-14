
from ..common cimport Table, Lookup, BoundaryMatrix, Index, Vector, Set

from libcpp.vector cimport vector
from libcpp.set cimport set
from libcpp.map cimport map

cdef extern from "LinBoxMethods.h":
	Vector LanczosKernelSample(Index coboundary, int M, int N, int p, int maxTries) except +
	Set ComputePercolationEvents(Table addition, Table multiplication, Lookup negation, Lookup inversion, BoundaryMatrix Boundary, Index breaks, int cellCount, int topDimension, int homology)
