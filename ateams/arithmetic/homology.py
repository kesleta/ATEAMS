
import numpy as np
import cython


@cython.boundscheck(False)
@cython.wraparound(False)
def reindexSparseBoundaryMatrix(
		cycles, boundary, flatBoundary, homology, tranches, reindexed, targetIndices,
		zeroedTargetIndices
	):
	filtration = np.arange(tranches[homology+1][1])

	# Evaluate the cochain on the complex to determine the (un)satisfied
	# plaquettes.
	satisfiedIndices = (cycles == 0).nonzero()[0]
	unsatisfiedIndices = (cycles > 0).nonzero()[0]

	# Shuffle the satisfied indices, then shuffle the indices in the filtration.
	low, high = tranches[homology]
	m = len(satisfiedIndices)

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
	target = reindexed[homology][filtered].tolist()

	_reindexer = np.array([filtered, zeroedTargetIndices])
	reindexer = _reindexer[:, _reindexer[0].argsort()][1]+low
	higher = reindexer[boundary[homology+1]].tolist()

	flattened = flatBoundary[:low] + target + higher

	return low, filtration, flattened, shuffledIndices, satisfiedIndices


@cython.boundscheck(False)
@cython.wraparound(False)
def computePersistencePairs(
		times, maxIndex, premarked, dimensions, maxDimension, tranches, rowSizes, fieldInverses, fieldCharacteristic,
		filtration, boundary
	):
	marked = set(premarked)
	dynamicCoefficients = { t: {} for t in range(maxIndex) }
	indexMap = { t: {} for t in range(maxIndex) }
	degrees = {}
	percolationEvents = set()

	for t in times:
		# All vertices initially have empty boundary, so they're all marked.
		cell = filtration[t]
		dim = dimensions[cell]
		low = tranches[dim-1][0]
		d = removePivotRows(boundary[cell], marked, dynamicCoefficients, indexMap, low, rowSizes[dim-1], fieldInverses, fieldCharacteristic)

		if not (d[0] > 0).any():
			marked |= {cell}
			degrees[cell] = t
		else:
			# Find the largest index in the linear combination (i.e the largest
			# degree, i.e. the "1 furthest to the right"), then drop all the
			# columns with coefficient 0.
			i = max(d[1])
			mask = (d[0] > 0)
			m = d[:, mask]
			dynamicCoefficients[i] = m
			indexMap[i] = dict(zip(m[1], m[0]))
			degrees[i] = t

	# Capture only the essential features; we don't care about anything else
	# (for now).
	for cell in marked:
		if dimensions[cell] != maxDimension-1: continue
		if not len(dynamicCoefficients[cell]): percolationEvents |= {degrees[cell]}

	return percolationEvents


@cython.boundscheck(False)
@cython.wraparound(False)
def removePivotRows(boundingChain, marked, T, indexMap, low, maxCells, inverses, p):
	# Create an extensible data structure for fast lookups and arithmetic. We're
	# trying out "big" arithmetic.
	_tossed = [b for b in boundingChain if b not in marked]
	tossed = np.array(_tossed)-low if len(_tossed) else []
	boundingChain = np.array(boundingChain)-low

	coefficients = np.zeros(maxCells, dtype=int)
	coefficients[boundingChain[::2]] = -1 % p
	coefficients[boundingChain[1::2]] = 1
	coefficients[tossed] = 0

	while (coefficients > 0).any():
		# Find the biggest coefficient and look up what's stored in the pivot
		# table. `i` is the index of the plaquette of greatest degree (added
		# latest) in the current boundary, `pivot` is the table of coefficients,
		# and `indexMap` is the mapping from degrees to coefficients in `T[i]`.
		# For example, if `i = 3` and `pivot = [[4, 2], [1, 3]]`, then `s = 1`
		# and we retrieve the entry in the first row of the second column as the
		# coefficient q.
		i = max((coefficients > 0).nonzero()[0])+low
		nonpivot = T[i]

		# If this column *is* a pivot, we're done.
		if not len(nonpivot): break

		# Otherwise, we're zeroing out the row. Given the compact information
		# stored, we simply fill in the values, do a subtraction, and we're out.
		s = indexMap[i]
		q = inverses[s[i]]

		pivotCoefficients = nonpivot[0]
		pivotIndices = nonpivot[1]-low
		coefficients[pivotIndices] = (coefficients[pivotIndices] - (q*pivotCoefficients)%p)%p

	gt = (coefficients > 0).nonzero()[0]
	# print(np.array([coefficients[gt], gt+low]))
	return np.array([coefficients[gt], gt+low])