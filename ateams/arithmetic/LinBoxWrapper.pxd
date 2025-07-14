
from ..common cimport INDEXFLAT, Index

cpdef Index LanczosKernelSample(INDEXFLAT coboundary, INDEXFLAT zeros, int faces, int columns, int p, int maxTries=*) except *
cpdef Index SubLanczosKernelSample(INDEXFLAT coboundary, INDEXFLAT zeroRows, INDEXFLAT zeroColumns, int p, int maxTries=*) except *
