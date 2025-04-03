
import numpy as np
from .common import FINT


def reindexSparseBoundaryMatrix(
		upperBoundary,
		flatBoundary,
		homology,
		tranches,
		reindexed,
		targetIndices,
		zeroedTargetIndices,
		cycles
	):
	filtration = np.arange(tranches[homology+1][1], dtype=FINT)

	# Evaluate the cochain on the complex to determine the (un)satisfied
	# plaquettes.
	satisfiedIndices = (cycles == 0).nonzero()[0]
	unsatisfiedIndices = (cycles > 0).nonzero()[0]

	# Shuffle the satisfied indices, then shuffle the indices in the filtration.
	low, high = tranches[homology]
	m = satisfiedIndices.shape[0]

	shuffledIndices = np.random.permutation(satisfiedIndices).astype(FINT)
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
