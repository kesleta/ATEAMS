
from libcpp.vector cimport vector

cdef extern from "LinBoxMethods.h":
	vector[int] LanczosKernelSample(vector[int] coboundary, int M, int N, int p, int maxTries);
