
import numpy as np
cimport numpy as np
import cython


@cython.boundscheck(False)
@cython.wraparound(False)
def reindexSparseBoundaryMatrix(
		np.ndarray[np.int_t, ndim=1] cycles,
		np.ndarray[np.int_t, ndim=2] upperBoundary,
		flatBoundary,
		Py_ssize_t homology,
		np.ndarray[np.int_t, ndim=2] tranches,
		np.ndarray[np.int_t, ndim=2] reindexed,
		np.ndarray[np.int_t, ndim=1] targetIndices,
		np.ndarray[np.int_t, ndim=1] zeroedTargetIndices
	):
	cdef np.ndarray[np.int_t, ndim=1] filtration;
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


# @cython.boundscheck(False)
# @cython.wraparound(False)
def computeGiantCyclePairs(
		np.ndarray[np.int_t, ndim=1] times,
		premarked,
		np.ndarray[np.int_t, ndim=1] dimensions,
		np.ndarray[np.int_t, ndim=1] fieldInverses,
		Py_ssize_t fieldCharacteristic,
		int maxDimension,
		int maxIndex,
		np.ndarray[np.int_t, ndim=1] filtration,
		boundary
	):
	marked = { *premarked }
	dynamicCoeffs = { t: {} for t in range(maxIndex) }
	cdef int O = len(filtration)
	cdef np.ndarray[np.int_t, ndim=1] degree = np.zeros(O, dtype=int)
	events = set()

	cdef Py_ssize_t _t;
	cdef Py_ssize_t t;
	cdef Py_ssize_t i;
	cdef Py_ssize_t cell;
	cdef int N = times.shape[0];

	for _t in range(N):
		t = times[_t]
		cell = filtration[t]
		chain = reducePivotRow(boundary[cell], marked, dynamicCoeffs, fieldInverses, fieldCharacteristic)

		if not len(chain):
			marked.add(cell)
			degree[cell] = t
		else:
			i = max(chain);
			dynamicCoeffs[i] = chain
			degree[i] = t
	
	cdef np.ndarray[np.int_t, ndim=1] remarked = np.array(list(marked));
	cdef int M = remarked.shape[0];
	cdef Py_ssize_t _s;

	for _s in range(M):
		cell = remarked[_s]
		if dimensions[cell] != maxDimension-1: continue
		if not len(dynamicCoeffs[cell]): events.add(degree[cell])
	
	return events


@cython.boundscheck(False)
@cython.wraparound(False)
def reducePivotRow(
		B,
		marked,
		dynamicCoeffs,
		np.ndarray[np.int_t, ndim=1] fieldInverses,
		int fieldCharacteristic
	):
	# Remove marked cells.
	includedCells = set(B) & marked

	cdef Py_ssize_t j;
	cdef int N = len(B);
	coefficients = {}

	for j in range(N):
		if B[j] in marked:
			coefficients[B[j]] = (-1)**(j+1) % fieldCharacteristic

	while len(includedCells):
		# Find the largest index over the included cells; this is our potential
		# pivot row.
		i = max(includedCells)
		maybe = dynamicCoeffs[i]

		# If `maybe` is a pivot row, then we're done, and we return the appropriate
		# coefficients.
		if not len(maybe): break

		# If `maybe` is *not* a pivot row, then we zero out the row. First, we
		# get the coefficient in `maybe` of the largest index (stored in `i`).
		# Then, we subtract off that multiple of the portion of the boundary
		# stored in `maybe` from `coefficients`.
		q = fieldInverses[int(maybe[i])-1]

		left, right = includedCells, set(maybe)

		_shared = {
			t: _smmod(coefficients[t], maybe[t], q, fieldCharacteristic)
			for t in left & right
		}

		_lonly = {
			t: coefficients[t] for t in left - right
		}

		_ronly = {
			t: _smmod(0, maybe[t], q, fieldCharacteristic) for t in right - left
		}

		coefficients = {} | _shared | _lonly | _ronly
		coefficients = {
			t: coefficients[t] for t in coefficients if coefficients[t] > 0
		}

		includedCells = set(coefficients)

	return coefficients


@cython.cfunc
cdef int _smmod(int a, int b, int c, int p):
	return _smod(a, _mmod(b, c, p), p)

@cython.cfunc
cdef int _smod(int a, int b, int p): return (a-b) % p

@cython.cfunc
cdef int _mmod(int a, int b, int p): return (a*b) % p
