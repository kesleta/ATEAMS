
from ..common cimport INDEXFLAT, Index, bool

cpdef Index ReducedKernelSample(INDEXFLAT coboundary, INDEXFLAT zeros, int faces, int columns, int p, bool verbose) except *
cpdef Index SubReducedKernelSample(INDEXFLAT coboundary, INDEXFLAT zeroRows, INDEXFLAT zeroColumns, int p, bool verbose) except *
