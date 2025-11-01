
from ..common cimport Index

cdef extern from "Sampling.h":
	Index ReducedKernelSample(Index coboundary, int M, int N, int p) except +
