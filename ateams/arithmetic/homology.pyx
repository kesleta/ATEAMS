
import numpy as np
cimport numpy as np
import cython

np.import_array()
ctypedef np.int64_t DTYPE_t
DTYPE = np.int64


@cython.boundscheck(False)
@cython.wraparound(False)
def reindexSparseBoundaryMatrix(
		np.ndarray[DTYPE_t, ndim=1] cycles,
		np.ndarray[DTYPE_t, ndim=2] upperBoundary,
		flatBoundary,
		Py_ssize_t homology,
		np.ndarray[DTYPE_t, ndim=2] tranches,
		np.ndarray[DTYPE_t, ndim=2] reindexed,
		np.ndarray[DTYPE_t, ndim=1] targetIndices,
		np.ndarray[DTYPE_t, ndim=1] zeroedTargetIndices
	):
	cdef np.ndarray[DTYPE_t, ndim=1] filtration;
	filtration = np.arange(tranches[homology+1][1])

	# Evaluate the cochain on the complex to determine the (un)satisfied
	# plaquettes.
	satisfiedIndices = (cycles == 0).nonzero()[0]
	unsatisfiedIndices = (cycles > 0).nonzero()[0]

	# Shuffle the satisfied indices, then shuffle the indices in the filtration.
	cdef Py_ssize_t low;
	cdef Py_ssize_t high;
	cdef Py_ssize_t m;
	low, high = tranches[homology]
	m = satisfiedIndices.shape[0]

	cdef np.ndarray[np.int_t, ndim=1] shuffledIndices;
	cdef np.ndarray[np.int_t, ndim=1] satisfiedCells;
	cdef np.ndarray[np.int_t, ndim=1] unsatisfiedCells;

	shuffledIndices = np.random.permutation(satisfiedIndices)
	satisfiedCells = targetIndices[shuffledIndices]
	unsatisfiedCells = targetIndices[unsatisfiedIndices]
	filtration[low:low+m] = satisfiedCells
	filtration[low+m:high] = unsatisfiedCells
	filtered = filtration[low:high]-low

	# Reconstruct a flattened boundary matrix. Shuffle the cells of target
	# dimension; afterward, we have to re-index the boundary matrix for the
	# cells of dimension one higher than the targets. This requires a mapping
	# from the previous index of each target cell to its new index. We can
	# find this mapping by creating a 2-d array with the shuffled (filtered)
	# indices as the top row and in-order zeroed-out indices as the bottom
	# row, then column-wise sorting this array by the values in the top row;
	# the second row of the sorted array is then the appropriate mapping.
	target = reindexed[filtered].tolist()

	_reindexer = np.array([filtered, zeroedTargetIndices])
	reindexer = _reindexer[:, _reindexer[0].argsort()][1]+low
	higher = reindexer[upperBoundary].tolist()

	flattened = flatBoundary[:low] + target + higher

	return low, filtration, flattened, shuffledIndices, satisfiedIndices


@cython.boundscheck(False)
@cython.wraparound(False)
def computeGiantCyclePairs(
		np.ndarray[DTYPE_t, ndim=1] times,
		premarked,
		np.ndarray[DTYPE_t, ndim=1] dimensions,
		np.ndarray[DTYPE_t, ndim=2] tranches,
		np.ndarray[DTYPE_t, ndim=1] fieldInverses,
		int fieldCharacteristic,
		int maxDimension,
		int maxIndex,
		np.ndarray[DTYPE_t, ndim=1] filtration,
		boundary,
		np.ndarray[DTYPE_t, ndim=1] degree
	):
	# Buckets for marked indices, dynamic coefficients (for storing 2d arrays
	# corresponding to chains), and indices (mapping cell degrees to indices in
	# `dynamicCoeffs`). The first data structure is required; the second data
	# structure tries to mitigate linear-time searches for max index coefficients.
	marked = { *premarked }
	dynamicCoeffs = { }
	dynamicIndices = { }
	dynamicSets = { }

	cdef Py_ssize_t _j;

	for _j in range(maxIndex):
		dynamicCoeffs[_j] = np.empty((2,0), dtype=DTYPE)
		dynamicIndices[_j] = {}
		dynamicSets[_j] = set();

	# Set for collecting the degrees of giant cycles.
	events = set()

	# We want fast loops, so this is how we do it.
	cdef Py_ssize_t _t;
	cdef Py_ssize_t t;
	cdef Py_ssize_t i;
	cdef Py_ssize_t cell;
	cdef Py_ssize_t dim;
	cdef int width;
	cdef np.ndarray[DTYPE_t, ndim=2] chain;
	cdef int N = times.shape[0];

	for _t in range(N):
		t = times[_t]
		cell = filtration[t]
		dim = dimensions[cell]
		width = tranches[dim-1,1]-tranches[dim-1,0] // 4
		
		chain = reducePivotRow(
			np.zeros((2, width), dtype=DTYPE),
			boundary[cell],
			marked,
			dynamicCoeffs,
			dynamicIndices,
			dynamicSets,
			fieldInverses,
			fieldCharacteristic
		)

		if chain.shape[1] < 1:
			marked.add(cell)
			degree[cell] = t
		else:
			i = _max(chain[1])
			dynamicCoeffs[i] = chain
			dynamicIndices[i] = dict(zip(chain[1], range(chain.shape[1])))
			dynamicSets[i] = set(chain[1])
			degree[i] = t
	
	cdef np.ndarray[np.int_t, ndim=1] remarked = np.array(list(marked));
	cdef int M = remarked.shape[0];
	cdef Py_ssize_t _s;

	for _s in range(M):
		cell = remarked[_s]
		if dimensions[cell] != maxDimension-1: continue
		if dynamicCoeffs[cell].shape[1] < 1: events.add(degree[cell])
	
	return events


@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray[DTYPE_t, ndim=2] reducePivotRow(
		np.ndarray[DTYPE_t, ndim=2] coefficients,
		boundary,
		marked,
		dynamicCoeffs,
		dynamicIndices,
		dynamicSets,
		np.ndarray[DTYPE_t, ndim=1] fieldInverses,
		int fieldCharacteristic
	):
	# Variable typing.
	cdef np.ndarray[DTYPE_t, ndim=1] nonzeroIndices;
	cdef Py_ssize_t _t, t, i, j, c, d, r, rindex, locator, a;
	cdef int L, R, S, q, _M, x, y, added;
	cdef int _N = len(boundary);

	# Fill the bucket of coefficients; we type this *outside* this function
	# because Cython doesn't like typing numpy arrays inline? Weird.
	for _t in range(_N):
		t = boundary[_t]
		if t not in marked: continue

		coefficients[0, _t] = (int)((-1)**(_t+1) % fieldCharacteristic)
		coefficients[1, _t] = t

	# Determine where the "occupied" and "free" indices are.
	cdef np.ndarray[DTYPE_t, ndim=1] occupied;
	cdef int _occupied = _countNonzero(coefficients)

	while _occupied > 0:
		# Find the largest index over the included cells; this is our potential
		# pivot row. Notably, we *cannot* include indices with coefficient zero.
		i = _maxNonzero(coefficients)
		nonpivot = dynamicCoeffs[i]

		# If the column *is* a pivot, then we've already reduced the row and we
		# return the coefficients.
		if nonpivot.shape[1] < 1: break

		# Otherwise, we determine the maximum index over those specified by the
		# chain and get its multiplicative inverse (over the finite field).
		# Because `dynamicIndices` is a Python `dict`, `s` has to remain untyped.
		# `s` is then just a dictionary mapping degrees to indices, so we use this
		# to access the coefficients.
		occupied = _nonzeroIndices(coefficients)
		coefficientIndices = dict()
		left = set()

		for d in range(occupied.shape[0]):
			r = occupied[d]
			coefficientIndices[coefficients[1,r]] = r
			left.add(coefficients[1,r])
		
		nonpivotIndices = dynamicIndices[i]
		locator = nonpivotIndices[i]
		q = fieldInverses[nonpivot[0,locator]]

		# Find the overlap in sets.
		right = set(nonpivot[1])
		shared = left & right
		rOnly = right - left
		L = len(left)
		R = len(rOnly)

		# Iterate over shared indices.
		for k in shared:
			lindex = coefficientIndices[k]
			rindex = nonpivotIndices[k]
			coefficients[0,lindex] = _smmod(
				coefficients[0,lindex],
				nonpivot[0,rindex],
				q,
				fieldCharacteristic
			)

		# Now, we can place values in any index that is zeroed out. We iterate
		# over `nonpivot` to check whether it's been added already.
		added = 0;
		a = 0;

		for a in range(coefficients.shape[1]):
			if added == R: break
			if coefficients[0,a] < 1:
				rindex = nonpivotIndices[rOnly.pop()]
				coefficients[1,a] = nonpivot[1, rindex]
				coefficients[0,a] = _smmod(
					0,
					nonpivot[0,rindex],
					q,
					fieldCharacteristic
				)

				added += 1
		
		_occupied = _countNonzero(coefficients)
	
	return _sliceMatrix(coefficients)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cfunc
cdef int _countNonzero(np.ndarray[DTYPE_t, ndim=2] A):
	cdef Py_ssize_t j = 0;
	cdef int N = A.shape[1];
	cdef int l = 0;

	for j in range(N):
		if A[0,j] > 0: l += 1

	return l


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cfunc
cdef np.ndarray[DTYPE_t, ndim=1] _nonzeroIndices(np.ndarray[DTYPE_t, ndim=2] A):
	cdef Py_ssize_t j = 0;
	cdef int N = A.shape[1];
	cdef int l = 0;
	cdef np.ndarray[DTYPE_t, ndim=1] indices = np.zeros(N, dtype=DTYPE);

	for j in range(N):
		if A[0,j] > 0:
			indices[l] = j;
			l += 1

	return indices[:l]


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cfunc
cdef np.ndarray[DTYPE_t, ndim=2] _sliceMatrix(np.ndarray[DTYPE_t, ndim=2] A):
	cdef Py_ssize_t j;
	cdef Py_ssize_t t;
	cdef np.ndarray[DTYPE_t, ndim=1] nonzero = _nonzeroIndices(A);
	cdef int N = nonzero.shape[0];
	cdef np.ndarray[DTYPE_t, ndim=2] resized = np.empty((2, N), dtype=DTYPE);

	for j in range(N):
		t = nonzero[j]
		resized[0,j] = A[0,t]
		resized[1,j] = A[1,t]

	return resized


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cfunc
cdef int _maxNonzero(np.ndarray[DTYPE_t, ndim=2] A):
	cdef Py_ssize_t j = 0;
	cdef int N = A.shape[1];
	cdef l = 0;

	for j in range(N):
		if A[0,j] > 0 and A[1,j] > l: l = A[1,j]

	return l


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cfunc
cdef int _max(np.ndarray[DTYPE_t, ndim=1] A):
	cdef Py_ssize_t j = 0;
	cdef int N = A.shape[0];
	cdef int l = A[0];

	for j in range(N):
		if A[j] > l: l = A[j]

	return l


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cfunc
cdef int _smmod(int a, int b, int c, int p):
	return _smod(a, _mmod(b, c, p), p)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cfunc
cdef int _smod(int a, int b, int p): return _mod(a-b, p)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cfunc
cdef int _mmod(int a, int b, int p): return _mod(a*b, p)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cfunc
@cython.cdivision(True)
cdef int _mod(int a, int b):
	cdef int r = a % b;

	if r < 0: return r + b;
	return r;
