
import numpy as np
from ..common import FINT


def reindexSparseBoundaryMatrix(
		homology,
		tranches,
		targetIndices,
		cycles
	):
	filtration = np.arange(tranches[homology+1][1])

	# Evaluate the cochain on the complex to determine the (un)satisfied
	# plaquettes.
	satisfiedIndices = (cycles == 0).nonzero()[0]
	unsatisfiedIndices = (cycles > 0).nonzero()[0]

	# Shuffle the satisfied indices, then shuffle the indices in the filtration.
	low, high = tranches[homology]
	m = satisfiedIndices.shape[0]

	shuffledIndices = np.random.permutation(satisfiedIndices)
	satisfiedCells = targetIndices[shuffledIndices]
	unsatisfiedCells = targetIndices[unsatisfiedIndices]
	filtration[low:low+m] = satisfiedCells
	filtration[low+m:high] = unsatisfiedCells

	return filtration, shuffledIndices, satisfiedIndices
