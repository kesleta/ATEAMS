

from ..common cimport INDEXFLAT

from libcpp.vector cimport vector as Vector
from libcpp.set cimport set as Set

cpdef Vector[int] LanczosKernelSample(INDEXFLAT coboundary, INDEXFLAT zeros, int faces, int columns, int p, int maxTries=*) except *
cpdef Vector[int] SubLanczosKernelSample(INDEXFLAT coboundary, INDEXFLAT zeroRows, INDEXFLAT zeroColumns, int p, int maxTries=*) except *
cpdef Set[int] ComputePercolationEvents(INDEXFLAT boundary, INDEXFLAT filtration, int homology, int p, INDEXFLAT breaks) noexcept
