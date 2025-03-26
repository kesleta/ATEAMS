
# cython: language_level=3str, initializedcheck=False, c_api_binop_methods=True, nonecheck=False, profile=True, cdivision=True, wraparound=False, boundscheck=False
# define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
# distutils: define_macros=CYTHON_TRACE_NOGIL=1
# distutils: language=c++

import numpy as np
cimport numpy as np
from .common cimport FFINT, FLAT, TABLE

from libcpp cimport bool
from .Sparse cimport Matrix


cdef void addRows(
		TABLE A,
		int i,
		int j,
		TABLE addition,
		TABLE multiplication,
		FFINT r
	) noexcept nogil:
	cdef int N, k;
	cdef FFINT mult;

	N = A.shape[1];

	for k in range(N):
		mult = multiplication[A[i,k],r]
		A[j,k] = addition[mult,A[j,k]]
		

cdef void swapRows(
		TABLE A,
		int i,
		int j
	) noexcept nogil:
	cdef int k, N;
	cdef FFINT t;

	N = A.shape[1];

	for k in range(N):
		t = A[i,k];
		A[i,k] = A[j,k];
		A[j,k] = t;


cdef void invertRow(
		TABLE A,
		int i,
		FFINT q,
		TABLE multiplication
	) noexcept nogil:
	cdef int k, N;
	
	N = A.shape[1]

	for k in range(N):
		A[i,k] = multiplication[q,A[i,k]]


cdef FFINT pivotRow(FLAT c, int pivots) noexcept nogil:
	cdef int N, i;

	N = c.shape[0];

	for i in range(N):
		if c[i] > 0 and i > pivots-1: return i
	
	return -1


cdef bool ContainsNonzero(FLAT b) noexcept nogil:
	cdef int i, N;
	N = b.shape[0];

	for i in range(N):
		if b[i] > 0: return True

	return False


cdef int MinZero(TABLE A) noexcept nogil:
	cdef int i, N;
	N = A.shape[0];

	for i in range(N):
		if not ContainsNonzero(A[i]): return i

	# If there are no nonzero rows, then we return -1.
	return -1


cdef TABLE RREF(
		TABLE A,
		int AUGMENT,
		FLAT P,
		TABLE addition,
		TABLE subtraction,
		FLAT negation,
		TABLE multiplication,
		FLAT inverses
	) noexcept nogil:
	cdef int M, N, i, j, k, pivot, pivots, ratio;
	cdef FFINT q;

	M = A.shape[0] # rows
	N = A.shape[1] # columns
	pivots = 0;

	# Compute REF.
	for i in range(AUGMENT):
		# Find the row with a pivot in this column (if one exists). If it doesn't,
		# continue; otherwise, mark and invert the row.
		pivot = pivotRow(A[:,i], pivots)
		if pivot < 0: continue

		# Store the pivot column and swap row entries.
		P[pivots] = i
		swapRows(A, pivots, pivot)

		# Invert the pivot row and increment the number of pivots.
		q = inverses[A[pivots,i]];
		invertRow(A, pivots, q, multiplication)
		pivots += 1;

		# Eliminate rows below.
		for k in range(pivots, M):
			ratio = negation[A[k,i]]
			addRows(A, pivots-1, k, addition, multiplication, ratio)

	# Compute RREF.
	cdef int l, m, r, column;
	l = pivots-1

	# Working backwards from the last nonzero row...
	while l > -1:
		# ... check that we're not getting an out-of-bounds error...
		if l > 0 and P[l] < 1:
			l -= 1;
			continue;

		# find the column containing the pivot element...
		column = P[l]

		# and, iterating over the rows *above* row `r`, eliminate the nonpivot
		# entries there...
		for m in range(l):
			ratio = negation[A[m, column]]
			addRows(A, l, m, addition, multiplication, ratio)

		# ... decrement the counter.
		l -= 1;
	
	return A


cpdef np.ndarray[FFINT, ndim=2] KernelBasis(
		TABLE coboundary,
		int AUGMENT,
		FLAT P,
		TABLE addition,
		TABLE subtraction,
		FLAT negation,
		TABLE multiplication,
		FLAT inverses
	):
	cdef TABLE inversion, reduced, superreduced;
	cdef Matrix M = Matrix(coboundary, addition, negation, multiplication, inverses)
	cdef int minzero;

	M.RREF(AUGMENT);
	minzero = M.HighestZeroRow(AUGMENT);
	inversion = M.ToArray();
	reduced = inversion[:,AUGMENT:];
	superreduced = reduced[minzero:]

	return np.asarray(superreduced)

