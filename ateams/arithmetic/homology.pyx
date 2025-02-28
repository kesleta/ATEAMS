
# cython: profile=True, boundscheck=False, wraparound=False, cdivision=True
# distutils: language=C

import numpy as np
cimport numpy as np
import cython

np.import_array()
ctypedef np.int64_t DTYPE_t
DTYPE = np.int64


cpdef computeGiantCyclePairs(
		DTYPE_t [:] times,
		premarked,
		DTYPE_t [:] dimensions,
		DTYPE_t [:,:] tranches,
		DTYPE_t [:] fieldInverses,
		DTYPE_t fieldCharacteristic,
		int maxDimension,
		int maxIndex,
		DTYPE_t [:,:] empty,
		DTYPE_t [:,:] chain,
		DTYPE_t [:] indices,
		DTYPE_t [:,:] addition,
		DTYPE_t [:,:] subtraction,
		DTYPE_t [:,:] multiplication,
		DTYPE_t [:] powers,
		DTYPE_t [:] filtration,
		boundary,
		DTYPE_t [:] degree
	):
	# Buckets for marked indices, dynamic coefficients (for storing 2d arrays
	# corresponding to chains), and indices (mapping cell degrees to indices in
	# `dynamicCoeffs`). The first data structure is required; the second data
	# structure tries to mitigate linear-time searches for max index coefficients.
	marked = set(premarked)
	dynamicCoeffs = { }
	dynamicIndices = { }
	dynamicSets = { }

	cdef int _j;

	for _j in range(maxIndex):
		dynamicCoeffs[_j] = empty
		dynamicIndices[_j] = {}
		dynamicSets[_j] = set();

	# Set for collecting the degrees of giant cycles.
	events = set()

	# We want fast loops, so this is how we do it.
	cdef int _t, t, i, cell;
	cdef int N = times.shape[0];
	cdef DTYPE_t [:,:] reduced;

	for _t in range(N):
		t = times[_t]
		cell = filtration[t]
		
		chain = reducePivotRow(
			chain,
			boundary[cell],
			marked,
			indices,
			dynamicCoeffs,
			dynamicIndices,
			dynamicSets,
			fieldInverses,
			fieldCharacteristic,
			addition,
			subtraction,
			multiplication,
			powers
		)

		reduced = _sliceMatrix(chain, indices)

		if reduced.shape[1] < 1:
			marked.add(cell)
			degree[cell] = t
		else:
			i = _max(reduced[1])
			dynamicCoeffs[i] = reduced
			dynamicIndices[i] = dict(zip(reduced[1], range(reduced.shape[1])))
			dynamicSets[i] = set(reduced[1])
			degree[i] = t

		chain = _zeroOut(chain)
	
	cdef DTYPE_t [:,:] unmarked;

	for cell in marked:
		if dimensions[cell] != maxDimension-1: continue

		unmarked = dynamicCoeffs[cell]
		if unmarked.shape[1] < 1: events.add(degree[cell])
	
	return events


cdef DTYPE_t [:,:] reducePivotRow(
		DTYPE_t [:,:] coefficients,
		boundary,
		marked,
		DTYPE_t [:] indices,
		dynamicCoeffs,
		dynamicIndices,
		dynamicSets,
		DTYPE_t [:] fieldInverses,
		DTYPE_t fieldCharacteristic,
		DTYPE_t [:,:] addition,
		DTYPE_t [:,:] subtraction,
		DTYPE_t [:,:] multiplication,
		DTYPE_t [:] powers
	) noexcept:
	# Variable typing.
	cdef int _t, t, i, j, c, d, r, lindex, rindex, locator, a
	cdef int L, R, S, q, _M, x, y, added
	cdef int _N = len(boundary)
	cdef DTYPE_t [:] occupied
	cdef DTYPE_t [:,:] nonpivot
	cdef DTYPE_t inverse = -1
	cdef DTYPE_t mul, sub

	# Fill the bucket of coefficients; we type this *outside* this function
	# because Cython doesn't like typing numpy arrays inline? Weird.
	for _t in range(_N):
		t = boundary[_t]
		if t not in marked: continue

		coefficients[0,_t] = powers[_t]
		coefficients[1,_t] = t

	# Determine where the "occupied" and "free" indices are.
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
		occupied = _nonzeroIndices(coefficients, indices)
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
			sub = coefficients[0,lindex]
			mul = multiplication[q,nonpivot[0,rindex]]
			coefficients[0,lindex] = subtraction[sub,mul]

		# Now, we can place values in any index that is zeroed out. We iterate
		# over `nonpivot` to check whether it's been added already.
		added = 0;
		a = 0;

		for a in range(coefficients.shape[1]):
			if added == R: break
			if coefficients[0,a] < 1:
				rindex = nonpivotIndices[rOnly.pop()]
				mul = multiplication[q,nonpivot[0,rindex]]
				
				coefficients[0,a] = subtraction[0,mul]
				coefficients[1,a] = nonpivot[1,rindex]

				added += 1

		_occupied = _countNonzero(coefficients)

	return coefficients


@cython.cfunc
cdef DTYPE_t [:,:] _zeroOut(DTYPE_t [:,:] A) noexcept:
	cdef int j, N;
	N = A.shape[1]

	for j in range(N):
		A[0,j] = 0
		A[1,j] = 0

	return A


@cython.cfunc
cdef int _countNonzero(DTYPE_t [:,:] A) noexcept:
	cdef int j = 0;
	cdef int N = A.shape[1];
	cdef int l = 0;
	cdef int tripleZeros = 0;

	for j in range(N):
		if A[0,j] > 0: l += 1

		# Don't want to traverse the whole array, that wastes a lot of time
		if j > 0:
			if A[1,j-1] == 0 and A[1,j] == 0:
				tripleZeros += 1

		if tripleZeros > 4: break

	return l


@cython.cfunc
cdef DTYPE_t [:] _nonzeroIndices(DTYPE_t [:,:] A, DTYPE_t [:] indices) noexcept:
	cdef int j = 0;
	cdef int N = A.shape[1];
	cdef int l = 0;
	cdef int tripleZeros = 0;

	for j in range(N):
		if A[0,j] > 0:
			indices[l] = j;
			l += 1

		# Don't want to traverse the whole array, that wastes a lot of time
		if j > 0:
			if A[1,j-1] == 0 and A[1,j] == 0:
				tripleZeros += 1

		if tripleZeros > 4: break

	return indices[:l]


@cython.cfunc
cdef DTYPE_t [:,:] _sliceMatrix(DTYPE_t [:,:] A, DTYPE_t [:] indices) noexcept:
	cdef int j, t;
	cdef DTYPE_t [:] nonzero = _nonzeroIndices(A, indices);
	cdef int N = nonzero.shape[0];
	cdef DTYPE_t [:,:] resized = np.empty((2, N), dtype=DTYPE);

	for j in range(N):
		t = nonzero[j]
		resized[0,j] = A[0,t]
		resized[1,j] = A[1,t]

	return resized


@cython.cfunc
cdef DTYPE_t _maxNonzero(DTYPE_t [:,:] A) noexcept:
	cdef int j = 0;
	cdef int N = A.shape[1];
	cdef DTYPE_t l = 0;
	cdef int tripleZeros = 0;

	for j in range(N):
		if A[0,j] > 0 and A[1,j] > l: l = A[1,j]

		# Don't want to traverse the whole array, that wastes a lot of time
		if j > 0:
			if A[1,j-1] == 0 and A[1,j] == 0:
				tripleZeros += 1

		if tripleZeros > 4: break

	return l




@cython.cfunc
cdef DTYPE_t _max(DTYPE_t [:] A) noexcept:
	cdef DTYPE_t j = 0;
	cdef DTYPE_t N = A.shape[0];
	cdef DTYPE_t l = A[0];

	for j in range(N):
		if A[j] > l: l = A[j]

	return l
