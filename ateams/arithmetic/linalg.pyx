
# cython: language_level=3str, initializedcheck=False, c_api_binop_methods=True, nonecheck=False, profile=True, cdivision=True, boundscheck=False, wraparound=False
# define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
# distutils: language=c++

import numpy as np
cimport numpy as np
from .common cimport FFINT, FLAT, TABLE
from .Sparse cimport Matrix


from cython.parallel import prange
from libcpp.vector cimport vector as Vector
from libcpp cimport bool
from libc.stdlib cimport rand, srand
from libc.time cimport time

# Seed the RNG.
srand(time(NULL));

cdef void addRows(
		TABLE A,
		int i,
		int j,
		int start,
		int stop,
		TABLE addition,
		TABLE multiplication,
		FFINT r
	) noexcept nogil:
	cdef int N, k;
	cdef FFINT mult;

	for k in range(start, stop):
		mult = multiplication[A[i,k],r]
		A[j,k] = addition[mult,A[j,k]]


cdef void parallelAddRows(
		TABLE A,
		int i,
		int j,
		int start,
		int stop,
		TABLE addition,
		TABLE multiplication,
		FFINT r,
		int cores,
		int minBlockSize,
		int maxBlockSize
	) noexcept nogil:
	cdef int block, blocks, MINCOL, MAXCOL, k, b, N, _q;
	cdef Vector[int] BLOCK;
	cdef Vector[Vector[int]] BLOCKS;
	
	# Compute the block size *of the current submatrix*, then schedule vertical
	# block computation appropriately: we have bounds on the "block size,"
	# so no thread is given too much (or too little) work to do. Choosing the
	# block sizes may be more of an art than a science, though.
	_q = (stop-start)//cores;
	block = min(max(minBlockSize, _q), min(maxBlockSize, _q))
	blocks = ((stop-start)//block)+1
	BLOCKS = Vector[Vector[int]]();

	for b in range(blocks):
		MINCOL = start+b*block;
		MAXCOL = start+(b+1)*block if start+(b+1)*block <= stop else stop
		
		if MINCOL >= stop: break

		BLOCK = Vector[int]();
		BLOCK.push_back(MINCOL);
		BLOCK.push_back(MAXCOL);
		BLOCKS.push_back(BLOCK);

	N = BLOCKS.size();

	# Compute block-wise.
	for k in prange(N, schedule="static", nogil=True, num_threads=cores, chunksize=1):
		MINCOL = BLOCKS[k][0]
		MAXCOL = BLOCKS[k][1]

		addRows(A, i, j, MINCOL, MAXCOL, addition, multiplication, r);


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
	
	N = A.shape[1];

	for k in range(N):
		A[i,k] = multiplication[q,A[i,k]]


cdef FFINT pivotRow(FLAT c, int pivots) noexcept nogil:
	cdef int N, i;

	N = c.shape[0];

	for i in range(N):
		if c[i] > 0 and i > pivots-1: return i
	
	return -1


cdef int ContainsNonzero(FLAT b) noexcept nogil:
	cdef int i, N;
	N = b.shape[0];

	for i in range(N):
		if b[i] > 0: return 1

	return 0


cdef int MinZero(TABLE A) noexcept nogil:
	cdef int i, N;
	N = A.shape[0];

	for i in range(N):
		if ContainsNonzero(A[i]) < 1: return i

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
		FLAT inverses,
		bool parallel,
		int minBlockSize,
		int maxBlockSize,
		int cores,
		str schedule
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

			if not parallel: addRows(A, pivots-1, k, i, N, addition, multiplication, ratio);
			else: parallelAddRows(A, pivots-1, k, i, N, addition, multiplication, ratio, cores, minBlockSize, maxBlockSize);

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

			if not parallel: addRows(A, l, m, column, N, addition, multiplication, ratio);
			else: parallelAddRows(A, l, m, column, N, addition, multiplication, ratio, cores, minBlockSize, maxBlockSize);

		# ... decrement the counter.
		l -= 1;
	
	return A


cpdef TABLE KernelBasis(
		FLAT P,
		TABLE addition,
		TABLE subtraction,
		FLAT negation,
		TABLE multiplication,
		FLAT inverses,
		TABLE coboundary,
		int AUGMENT,
		bool parallel,
		int minBlockSize,
		int maxBlockSize,
		int cores,
		str schedule,
	):
	cdef TABLE inversion, reduced, superreduced;
	cdef int minzero;

	inversion = RREF(coboundary, AUGMENT, P, addition, subtraction, negation, multiplication, inverses, parallel, minBlockSize, maxBlockSize, cores, schedule)
	minzero = MinZero(inversion)
	reduced = inversion[:,AUGMENT:];
	superreduced = reduced[minzero:]

	return superreduced


cdef int Random(int MAX) except -1:
	return (rand() % MAX);


cdef FLAT LinearCombination(
		TABLE basis,
		FLAT coefficients,
		FLAT store,
		TABLE addition,
		TABLE multiplication,
		int FIELD,
		FLAT result,
		bool parallel
	) noexcept:
	cdef int i, j, N, M;
	cdef FFINT mult;

	M = basis.shape[0];
	N = basis.shape[1];

	# Pick (uniform random) coefficients in the field.
	for i in range(M):
		store[i] = coefficients[Random(FIELD)]

	for j in range(N):
		for i in range(M):
			mult = multiplication[store[i], basis[i,j]]
			result[j] = addition[result[j], mult]

	return result


cpdef np.ndarray[FFINT, ndim=1, negative_indices=False, mode="c"] SampleFromKernel(
		int FIELD,
		FLAT PIVOTS,
		FLAT store,
		FLAT result,
		TABLE addition,
		TABLE subtraction,
		FLAT negation,
		TABLE multiplication,
		FLAT inverses,
		bool parallel,
		int minBlockSize,
		int maxBlockSize,
		int cores,
		str schedule,
		TABLE coboundary,
		int AUGMENT
	):
	cdef TABLE basis = KernelBasis(PIVOTS, addition, subtraction, negation, multiplication, inverses, coboundary, AUGMENT, parallel, minBlockSize, maxBlockSize, cores, schedule);
	return np.asarray(LinearCombination(basis, negation, store, addition, multiplication, FIELD, result, parallel))



cpdef TABLE SparseKernelBasis(
		FLAT P,
		TABLE addition,
		TABLE subtraction,
		FLAT negation,
		TABLE multiplication,
		FLAT inverses,
		TABLE coboundary,
		int AUGMENT,
		bool parallel,
		int minBlockSize,
		int maxBlockSize,
		int cores,
		str schedule
	):
	cdef TABLE inversion, reduced, superreduced;
	cdef Matrix M = Matrix(coboundary, addition, negation, multiplication, inverses, parallel, minBlockSize, maxBlockSize, cores, schedule)
	cdef int minzero;

	M.RREF(AUGMENT);
	minzero = M.HighestZeroRow(AUGMENT);
	inversion = M.ToArray();
	reduced = inversion[:,AUGMENT:];
	superreduced = reduced[minzero:]

	return superreduced


cpdef np.ndarray[FFINT, ndim=1, negative_indices=False, mode="c"] SparseSampleFromKernel(
		int FIELD,
		FLAT PIVOTS,
		FLAT store,
		FLAT result,
		TABLE addition,
		TABLE subtraction,
		FLAT negation,
		TABLE multiplication,
		FLAT inverses,
		bool parallel,
		int minBlockSize,
		int maxBlockSize,
		int cores,
		str schedule,
		TABLE coboundary,
		int AUGMENT
	):
	cdef TABLE basis = SparseKernelBasis(PIVOTS, addition, subtraction, negation, multiplication, inverses, coboundary, AUGMENT, parallel, minBlockSize, maxBlockSize, cores, schedule);
	return np.asarray(LinearCombination(basis, negation, store, addition, multiplication, FIELD, result, parallel))


