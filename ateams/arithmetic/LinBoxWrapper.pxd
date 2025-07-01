

from ..common cimport FFINT, TABLE, MINT, INDEXFLAT
from .LinBoxMethods cimport LanczosKernelSample as _LanczosKernelSample

from libcpp.vector cimport vector as Vector

cpdef Vector[int] LanczosKernelSample(INDEXFLAT coboundary, INDEXFLAT zeros, int faces, int columns, int field, int maxTries=*) noexcept

