
# distutils: language=c++

from ..common cimport INDEXFLAT, PersistencePairs

cpdef PersistencePairs PHATComputePersistencePairs(INDEXFLAT boundary, INDEXFLAT filtration, int homology, INDEXFLAT breaks) noexcept