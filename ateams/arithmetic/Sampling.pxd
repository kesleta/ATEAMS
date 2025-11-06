
from ..common cimport Index, bool

cdef extern from "Sampling.h":
	Index ReducedKernelSample(Index coboundary, int M, int N, int p, bool verbose) except +
