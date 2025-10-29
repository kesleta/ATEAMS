
from ..common cimport INDEXFLAT, Index

cpdef Index ReducedKernelSample(INDEXFLAT coboundary, INDEXFLAT zeros, int faces, int columns, int p) except *
cpdef Index SubReducedKernelSample(INDEXFLAT coboundary, INDEXFLAT zeroRows, INDEXFLAT zeroColumns, int p) except *
