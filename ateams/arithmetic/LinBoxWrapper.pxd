

from ..common cimport FFINT, TABLE, MINT, INDEXFLAT
from .LinBoxMethods cimport LanczosKernelSample as _LanczosKernelSample

from libcpp.vector cimport vector as Vector
from libcpp.set cimport set as Set

cpdef Vector[int] LanczosKernelSample(INDEXFLAT coboundary, INDEXFLAT zeros, int faces, int columns, int field, int maxTries=*) noexcept
cpdef Set[int] ComputePercolationEvents(INDEXFLAT boundary, INDEXFLAT filtration, int homology, int p, INDEXFLAT breaks) noexcept
