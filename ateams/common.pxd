
# distutils: language=c++

from libcpp.set cimport set as AnySet
from libcpp.vector cimport vector as AnyVector
from libcpp.map cimport map as AnyMap
from libcpp cimport bool

import numpy as np
cimport numpy as np

# Cython memoryview types. We should be switching to chars (or, really, abandoning
# these types altogether).
ctypedef np.int64_t MINT
ctypedef np.int16_t FFINT
ctypedef FFINT[:] FLAT
ctypedef FFINT[:,:] TABLE
ctypedef MINT[:,:] INDEXTABLE
ctypedef MINT[:] INDEX
ctypedef MINT[::1] INDEXFLAT
ctypedef FFINT[::1] FLATCONTIG
ctypedef FFINT[:,::1] TABLECONTIG

ctypedef AnyVector[int] Vector;
ctypedef AnyVector[int] Index;
ctypedef AnyVector[char] Lookup;
ctypedef AnyVector[AnyVector[char]] Table;
ctypedef AnyMap[int, char] Column;
ctypedef AnyVector[Column] BoundaryMatrix;
ctypedef AnyVector[Vector] Pairs;
ctypedef AnySet[int] Set;


cdef inline Vector Vectorize(INDEXFLAT A) noexcept:
	cdef int i, L;
	cdef Vector B;

	L = A.shape[0];
	B = Vector(L);

	for i in range(L): B[i] = A[i];

	return B;

