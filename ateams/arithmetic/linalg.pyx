
# cython: language_level=3str, initializedcheck=False, c_api_binop_methods=True, nonecheck=False, profile=True, cdivision=True
# cython: boundscheck=False, wraparound=False
# cython: binding=True, linetrace=True
# define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
# distutils: language=c++

import numpy as np
cimport numpy as np
from .common cimport FFINT, FLAT, TABLE
from .Sparse cimport Matrix


from cython.parallel cimport prange
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
	) noexcept:
	cdef int k;
	cdef FFINT mult;

	for k in range(start, stop):
		mult = multiplication[A[i,k],r]
		A[j,k] = addition[mult,A[j,k]]


cdef void ParallelAddRows(
		TABLE A,
		int i,
		int j,
		int start,
		int stop,
		TABLE addition,
		TABLE multiplication,
		FFINT r
	) noexcept nogil:
	cdef int k;
	cdef FFINT mult;

	for k in range(start, stop):
		mult = multiplication[A[i,k],r]
		A[j,k] = addition[mult,A[j,k]]


cdef Vector[Vector[int]] computeBlockSchema(
		int start,
		int stop,
		int cores,
		int minBlockSize,
		int maxBlockSize
	) noexcept:
	cdef int block, blocks, MINCOL, MAXCOL, b, _q;
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
	
	return BLOCKS


cdef void swapRows(TABLE A, int i, int j) noexcept:
	cdef int k, N;
	cdef FFINT t;

	N = A.shape[1];

	for k in range(N):
		t = A[i,k];
		A[i,k] = A[j,k];
		A[j,k] = t;


cdef void invertRow(TABLE A, int i, FFINT q, TABLE multiplication) noexcept:
	cdef int k, N;
	
	N = A.shape[1];

	for k in range(N):
		A[i,k] = multiplication[q,A[i,k]]


cdef FFINT pivotRow(FLAT c, int pivots) noexcept:
	cdef int N, i;

	N = c.shape[0];

	for i in range(N):
		if c[i] > 0 and i > pivots-1: return i
	
	return -1


cdef int ContainsNonzero(FLAT b, int AUGMENT) noexcept:
	cdef int i, N;

	if AUGMENT < 0: N = b.shape[0];
	else: N = AUGMENT

	for i in range(N):
		if b[i] > 0: return 1

	return 0


cdef int MinZero(TABLE A, int AUGMENT) noexcept:
	cdef int i, N;
	N = A.shape[0];

	for i in range(N):
		if ContainsNonzero(A[i], AUGMENT) < 1: return i

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
		int cores
	) noexcept:
	cdef int _t, t, M, N, B, i, k, row, block, pivot, pivots, ratio, start, stop;
	cdef FFINT q;
	cdef Vector[Vector[int]] BLOCKS;
	cdef Vector[int] negations = Vector[int](A.shape[0])

	# Catch for full-matrix reduction.
	if AUGMENT < 0: AUGMENT = A.shape[1];

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

		# Use parallelism here: we can eliminate chunks of the matrix at a time.
		if not parallel:
			for k in range(pivots, M):
				ratio = negation[A[k,i]]
				addRows(A, pivots-1, k, i, N, addition, multiplication, ratio);
		else:
			# Compute the block schema (i.e. the contiguous sets of columns
			# assigned to each thread at this step), read the negations of each
			BLOCKS = computeBlockSchema(i, N, cores, minBlockSize, maxBlockSize);
			for _t in range(pivots, M): negations[_t] = negation[A[_t, i]]

			for block in prange(BLOCKS.size(), nogil=True, schedule="dynamic", num_threads=cores):
				start = BLOCKS[block][0];
				stop = BLOCKS[block][1];
				t = pivots;

				while t < M:
					ratio = negations[t];
					ParallelAddRows(A, pivots-1, t, start, stop, addition, multiplication, ratio)
					t = t+1
			
		

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
		
		# ... and, iterating over the rows *above* row `l`, eliminate the
		# nonzero entries there...
		if not parallel:
			for m in range(l):
				ratio = negation[A[m,column]]
				addRows(A, l, m, column, N, addition, multiplication, ratio);
		else:
			# Compute the block schema (i.e. the contiguous sets of columns
			# assigned to each thread at this step), read the negations of each
			BLOCKS = computeBlockSchema(column, N, cores, minBlockSize, maxBlockSize);
			for _t in range(l): negations[_t] = negation[A[_t, column]]

			for block in prange(BLOCKS.size(), nogil=True, schedule="dynamic", num_threads=cores):
				start = BLOCKS[block][0];
				stop = BLOCKS[block][1];
				m = 0;

				while m < l:
					ratio = negations[m]
					ParallelAddRows(A, l, m, start, stop, addition, multiplication, ratio)
					m = m+1

		# ... decrement the counter.
		l -= 1;
	
	return A


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


cdef TABLE _KernelBasis(
		FLAT P,
		TABLE empty,
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
		int cores
	):
	cdef TABLE inversion, reduced, superreduced, eliminated;
	cdef int minzero;

	inversion = RREF(coboundary, AUGMENT, P, addition, subtraction, negation, multiplication, inverses, parallel, minBlockSize, maxBlockSize, cores)
	minzero = MinZero(inversion, AUGMENT)

	# If there's trivial kernel (i.e. we're full-rank), return nothing. Otherwise,
	# RREF the kernel basis (or just send the kernel back as-is).
	if minzero < 0:
		eliminated = empty;
	else:
		minzero = MinZero(inversion, AUGMENT)
		reduced = inversion[:,AUGMENT:];
		superreduced = reduced[minzero:]
		eliminated = superreduced

	return eliminated


cpdef TABLE KernelBasis(
		FLAT P,
		TABLE empty,
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
		int cores
	):
	cdef TABLE inversion, reduced, superreduced, eliminated;
	cdef int minzero;

	inversion = RREF(coboundary, AUGMENT, P, addition, subtraction, negation, multiplication, inverses, parallel, minBlockSize, maxBlockSize, cores)
	minzero = MinZero(inversion, AUGMENT)

	# If there's trivial kernel (i.e. we're full-rank), return nothing. Otherwise,
	# RREF the kernel basis (or just send the kernel back as-is).
	if minzero < 0:
		eliminated = empty;
	else:
		minzero = MinZero(inversion, AUGMENT)
		reduced = inversion[:,AUGMENT:];
		superreduced = reduced[minzero:]
		eliminated = RREF(superreduced, -1, P, addition, subtraction, negation, multiplication, inverses, parallel, minBlockSize, maxBlockSize, cores)

	return eliminated


cpdef np.ndarray[FFINT, ndim=1, negative_indices=False, mode="c"] SampleFromKernel(
		int FIELD,
		FLAT PIVOTS,
		TABLE empty,
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
		TABLE coboundary,
		int AUGMENT
	):
	cdef TABLE basis = _KernelBasis(PIVOTS, empty, addition, subtraction, negation, multiplication, inverses, coboundary, AUGMENT, parallel, minBlockSize, maxBlockSize, cores);
	return np.asarray(LinearCombination(basis, negation, store, addition, multiplication, FIELD, result, parallel))


cdef TABLE _SparseKernelBasis(
		FLAT P,
		TABLE empty,
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
		int cores
	):
	cdef TABLE inversion, reduced, superreduced;
	cdef Matrix M = Matrix(coboundary, addition, negation, multiplication, inverses, parallel, minBlockSize, maxBlockSize, cores)
	cdef int minzero;

	M.RREF(AUGMENT);
	minzero = M.HighestZeroRow(AUGMENT);

	if minzero < 0:
		eliminated = empty;
	else:
		inversion = M.ToArray();
		reduced = inversion[:,AUGMENT:];
		superreduced = reduced[minzero:]
		eliminated = superreduced

	return eliminated


cpdef TABLE SparseKernelBasis(
		FLAT P,
		TABLE empty,
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
		int cores
	):
	cdef TABLE inversion, reduced, superreduced;
	cdef Matrix G, M = Matrix(coboundary, addition, negation, multiplication, inverses, parallel, minBlockSize, maxBlockSize, cores)
	cdef int minzero;

	M.RREF(AUGMENT);
	minzero = M.HighestZeroRow(AUGMENT);

	if minzero < 0:
		eliminated = empty;
	else:
		inversion = M.ToArray();
		reduced = inversion[:,AUGMENT:];
		superreduced = reduced[minzero:]
		G = Matrix(superreduced, addition, negation, multiplication, inverses, parallel, minBlockSize, maxBlockSize, cores)
		G.RREF()
		eliminated = G.ToArray()

	return eliminated


cpdef np.ndarray[FFINT, ndim=1, negative_indices=False, mode="c"] SparseSampleFromKernel(
		int FIELD,
		FLAT PIVOTS,
		TABLE empty,
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
		TABLE coboundary,
		int AUGMENT
	):
	cdef TABLE basis = _SparseKernelBasis(PIVOTS, empty, addition, subtraction, negation, multiplication, inverses, coboundary, AUGMENT, parallel, minBlockSize, maxBlockSize, cores);
	return np.asarray(LinearCombination(basis, negation, store, addition, multiplication, FIELD, result, parallel))


