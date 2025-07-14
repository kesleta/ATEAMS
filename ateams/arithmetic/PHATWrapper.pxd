
# distutils: language=c++

from ..common cimport INDEXFLAT, Vectorize, Pairs, Vector

cpdef Pairs PHATComputePersistencePairs(INDEXFLAT boundary, INDEXFLAT filtration, int homology, INDEXFLAT breaks) noexcept