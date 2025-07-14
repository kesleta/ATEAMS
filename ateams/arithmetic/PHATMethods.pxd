
from ..common cimport Index, Pairs

cdef extern from "PHATMethods.h":
	Pairs PHATComputePersistencePairs(Index boundary, Index filtration, int homology, Index breaks) noexcept

