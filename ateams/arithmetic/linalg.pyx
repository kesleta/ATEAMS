
# cython: language_level=3str, initializedcheck=False, c_api_binop_methods=True, nonecheck=False, profile=True, cdivision=True, wraparound=False, boundscheck=False
# define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
# distutils: language=c

import numpy as np
cimport numpy as np

ctypedef np.int64_t DTYPE_t
DTYPE = np.int64


cdef DTYPE_t[:] Tare(DTYPE_t[:] b) noexcept nogil:
	cdef int i, N;
	N = b.shape[0];

	for i in range(N): b[i] = 0;

	return b;


cdef DTYPE_t[:,:] RREF(
		DTYPE_t[:,:] A,
		DTYPE_t[:,:] I,
		DTYPE_t[:] P,
		DTYPE_t[:,:] addition,
		DTYPE_t[:,:] subtraction,
		DTYPE_t[:] negation,
		DTYPE_t[:,:] multiplication,
		DTYPE_t[:] inverses
	) noexcept nogil:
	cdef int M, N, i, j, k, pivot, pivots, ratio;
	cdef DTYPE_t q;

	M = A.shape[0] # rows
	N = A.shape[1] # columns
	pivots = 0;

	# Compute REF.
	for i in range(N):
		# Find the row with a pivot in this column (if one exists). If it doesn't,
		# continue; otherwise, mark and invert the row.
		pivot = pivotRow(A[:,i], pivots)
		if pivot < 0: continue

		# Store the pivot column and swap row entries.
		P[pivots] = i
		A = swapRows(A, pivots, pivot)

		# Invert the pivot row and increment the number of pivots.
		q = inverses[A[pivots,i]];
		A = invertRow(A, pivots, q, multiplication)
		pivots += 1;

		# Eliminate rows below.
		for k in range(pivots, M):
			ratio = negation[A[k,i]]
			A = addRows(A, pivots-1, k, addition, multiplication, ratio)

	# Compute RREF.
	cdef int l, m, r, column;

	# Working backwards from the last nonzero row...
	for l in range(pivots):
		# find the column containing the pivot element...
		r = pivots-l
		column = P[r]

		# and, iterating over the rows *above* row `r`, eliminate the nonpivot
		# entries there.
		for m in range(r):
			ratio = negation[A[m, column]]
			A = addRows(A, r, m, addition, multiplication, ratio)
		
	return A


cdef DTYPE_t[:,:] addRows(
		DTYPE_t[:,:] A,
		int i,
		int j,
		DTYPE_t[:,:] addition,
		DTYPE_t[:,:] multiplication,
		DTYPE_t r
	) noexcept nogil:
	cdef int N, k;
	cdef DTYPE_t mult;

	N = A.shape[1];

	for k in range(N):
		mult = multiplication[A[i,k],r]
		A[j,k] = addition[mult,A[j,k]]

	return A
		

cdef DTYPE_t[:,:] swapRows(
		DTYPE_t[:,:] A,
		int i,
		int j
	) noexcept nogil:
	cdef int k, N;
	cdef DTYPE_t t;

	N = A.shape[1];

	for k in range(N):
		t = A[i,k];
		A[i,k] = A[j,k];
		A[j,k] = t;

	return A


cdef DTYPE_t[:,:] invertRow(
		DTYPE_t[:,:] A,
		int i,
		DTYPE_t q,
		DTYPE_t[:,:] multiplication
	) noexcept nogil:
	cdef int k, N;
	
	N = A.shape[1]

	for k in range(N):
		A[i,k] = multiplication[q,A[i,k]]

	return A


cdef DTYPE_t pivotRow(DTYPE_t[:] c, int pivots) noexcept nogil:
	cdef int N, i;

	N = c.shape[0];

	for i in range(N):
		if c[i] > 0 and i > pivots-1: return i
	
	return -1


cpdef nullspace(
		DTYPE_t[:,:] coboundary,
		DTYPE_t[:,:] I,
		DTYPE_t[:] P,
		DTYPE_t[:,:] addition,
		DTYPE_t[:,:] subtraction,
		DTYPE_t[:] negation,
		DTYPE_t[:,:] multiplication,
		DTYPE_t[:] inverses
	):
	return np.asarray(RREF(coboundary, I, P, addition, subtraction, negation, multiplication, inverses))
