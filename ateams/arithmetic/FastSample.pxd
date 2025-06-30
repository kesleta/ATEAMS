
from libcpp.vector cimport vector

cdef extern from "FastSample.h":
	vector[int] FastSample(vector[vector[int]], int p);
