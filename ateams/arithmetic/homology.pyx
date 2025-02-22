
import numpy as np
import cython

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
		times, maxIndex, premarked, dimensions, maxDimension, fieldInverses, fieldCharacteristic,
		filtration, boundary
	):
	marked = set(premarked)
	dynamicCoefficients = { t: {} for t in range(maxIndex) }
	degrees = {}
	percolationEvents = set()

	for t in times:
		# All vertices initially have empty boundary, so they're all marked.
		cell = filtration[t]
		d = removePivotRows(boundary[cell], marked, dynamicCoefficients, fieldInverses, fieldCharacteristic)

		if not len(d):
			marked |= {cell}
			degrees[cell] = t
		else:
			i = max(d)
			dynamicCoefficients[i] = d
			degrees[i] = t

	# Capture only the essential features; we don't care about anything else
	# (for now).
	for cell in marked:
		if dimensions[cell] != maxDimension-1: continue
		if not len(dynamicCoefficients[cell]): percolationEvents |= {degrees[cell]}

	return percolationEvents

@cython.boundscheck(False)
@cython.wraparound(False)
def removePivotRows(boundingChain, marked, T, inverses, p):
	# Determine coefficients and remove unmarked ones.
	coefficients = { a: ((-1)**(j+1) % p) for j, a in enumerate(boundingChain) if a in marked }

	while len(coefficients):
		i = max(coefficients)
		pivot = T[i]

		# If this column *is* a pivot, we're done.
		if not len(pivot): break

		# Otherwise, we're zeroing out the row.
		q = inverses[pivot[i]]

		left, right = set(coefficients), set(pivot)
		shared = left & right

		Shared = {
			t: (coefficients[t] - (q*pivot[t])%p)%p for t in shared
		}

		Left = {
			t: coefficients[t] for t in left-right
		}

		Right = { t: (-(q*pivot[t])%p)%p for t in right-left }

		coefficients = {} | Shared | Left | Right
		coefficients = { t: coefficients[t] for t in coefficients if coefficients[t] > 0 }

	return coefficients