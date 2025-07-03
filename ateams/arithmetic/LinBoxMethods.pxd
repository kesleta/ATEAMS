
from libcpp.vector cimport vector
from libcpp.set cimport set

cdef extern from "LinBoxMethods.h":
	vector[int] LanczosKernelSample(vector[int] coboundary, int M, int N, int p, int maxTries);
	set[int] ComputePercolationEvents(vector[int] boundary, vector[int] filtration, int homology, int p, vector[int] breaks);
