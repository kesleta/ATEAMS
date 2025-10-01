
import numpy as np
from math import comb

from ..arithmetic import ComputePersistencePairs
from ..common import Matrices, FINT


class Bernoulli():
	_name = "Bernoulli"
	
	def __init__(self, C, dimension=1, **kwargs):
		"""
		Initializes classic Bernoulli percolation on the provided complex,
		detecting percolation in the `dimension`-1th homology group.

		Args:
			C (Complex): The `Complex` object on which we'll be running experiments.
			dimension (int=1): The dimension of cells on which we're percolating.
		"""
		# Object access. Set a field.
		self.complex = C
		self.dimension = dimension
		self._returns = 1
		self.field = 2

		# Set a phantom spins attribute so we don't break the Chain
		self.spins = None


		# Force-recompute the matrices for a different dimension; creates
		# a set of orientations for fast elementwise products.
		self.matrices = Matrices()
		self.matrices.full = self.complex.matrices.full

		boundary, coboundary = self.complex.recomputeBoundaryMatrices(dimension)
		self.matrices.boundary = boundary
		self.matrices.coboundary = coboundary

		# Useful values to have later.
		self.cellCount = len(self.complex.flattened)
		self.cells = len(self.complex.Boundary[self.dimension])
		self.faces = len(self.complex.Boundary[self.dimension-1])
		self.target = np.arange(self.complex.breaks[self.dimension], self.complex.breaks[self.dimension+1])

		# Premake the "occupied cells" array; change the dimension of the complex
		# to correspond to the provided dimension.
		self.rank = comb(len(self.complex.corners), self.dimension)
		self.nullity = len(self.complex.Boundary[self.dimension])

		# Delegates computation for persistence.
		self._delegateComputation()

		# Seed the random number generator.
		self.RNG = np.random.default_rng()

	
	def _delegateComputation(self):
		low, high = self.complex.breaks[self.dimension], self.complex.breaks[self.dimension+1]
		times = set(range(self.cellCount))

		def whittle(pairs):
			_births, _deaths = zip(*pairs)
			births = set(_births)
			deaths = set(_deaths)

			return set(
				e for e in times-(births|deaths)
				if low <= e < high
			)

		def persist(filtration):
			essential = ComputePersistencePairs(
				self.matrices.full, filtration, self.dimension, self.complex.breaks
			)

			return whittle(essential)

		self.persist = persist


	def _filtrate(self):
		"""
		Constructs a filtration based on the evaluation of the cochain.
		"""
		# Find which cubes get zeroed out (i.e. are sent to zero by the cocycle).
		low = self.complex.breaks[self.dimension]
		high = self.complex.breaks[self.dimension+1]

		filtration = np.arange(self.cellCount)
		shuffled = np.random.permutation(self.target)
		filtration[low:high] = shuffled

		return filtration, shuffled-low


	def _initial(self):
		"""
		Computes an initial state for the model's Complex.

		Returns:
			A numpy `np.array` representing a vector of spin assignments.
		"""
		return self.RNG.integers(
			0, high=self.field, dtype=FINT, size=self.faces
		)
	

	def _proposal(self, time):
		"""
		Proposal scheme for generalized invaded-cluster evolution on the
		random-cluster model.

		Args:
			time (int): Step in the chain.

		Returns:
			A numpy array representing a vector of spin assignments.
		"""
		# Construct the filtration and find the essential cycles.
		filtration, shuffled = self._filtrate()
		essential = self.persist(filtration)

		j = 0
		low = self.complex.breaks[self.dimension]

		occupied = np.zeros((self.rank, self.nullity))

		for t in essential:
			occupiedIndices = shuffled[:t-low]
			occupied[j,occupiedIndices] = 1
			j += 1

		return occupied
	

	def _assign(self, cocycle): pass

