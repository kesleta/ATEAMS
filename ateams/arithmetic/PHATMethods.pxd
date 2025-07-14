
from ..common cimport Index, PersistencePairs

cdef extern from "PHATMethods.h":
	PersistencePairs PHATComputePersistencePairs(Index boundary, Index filtration, int homology, Index breaks) noexcept

