
from libcpp.vector cimport vector as Vector

cdef extern from "PHATMethods.h":
	Vector[Vector[int]] PHATComputePersistencePairs(Vector[int] boundary, Vector[int] filtration, int homology, Vector[int] breaks) noexcept

