
import numpy as np
import phat
import warnings

from math import comb

from ..arithmetic import (
	ComputePercolationEvents, LanczosKernelSample, MatrixReduction,
	Persistence, KernelSample, FINT
)
from ..common import Matrices, TooSmallWarning
from .Model import Model


class Bernoulli(Model):
	name = "Bernoulli"
	
	def __init__(
			self, C, dimension=1, initial=None, stop=lambda: 1, LinBox=True, sparse=True,
			parallel=False, minBlockSize=32, maxBlockSize=64, cores=4
		):
		"""
		Initializes classic Bernoulli percolation on the provided complex,
		detecting percolation in the `dimension`-1th homology group.

		Args:
			C (Complex): The `Complex` object on which we'll be running experiments.
			dimension (int=1): The dimension of cells on which we're percolating.
			initial (np.array): A vector of spin assignments to components.
			stop (function): A function that returns the number of essential cycles
				found before sampling the next configuration.
			LinBox (bool=True): Uses fast LinBox routines instead of slow inbuilt
				ones. WARNING: using inbuilt methods may dramatically increase
				computation time.
			sparse (boolean): Should matrices be formatted sparsely? (Uses C/C++).
			parallel (boolean): Should matrix computations be done in parallel? (Uses C/C++).
			minBlockSize (int=32): If `parallel` is truthy, this is the smallest
				number of columns processed in parallel.
			maxBlockSize (int=64): If `parallel` is truthy, this is the largest
				number of columns processed in parallel.
			cores (int=4): Number of available CPUs/cores/threads on the machine.
		"""
		# Object access.
		self.complex = C
		self.dimension = dimension
		self.stop = stop


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


		# Delegates computation for persistence and cocycle sampling.
		self._delegateComputation(LinBox, sparse, parallel, minBlockSize, maxBlockSize, cores)


		# Seed the random number generator.
		self.RNG = np.random.default_rng()


		# If no initial spin configuration is passed, initialize.
		if not initial: self.spins = self.initial()
		else: self.spins = (initial%self.complex.field).astype(FINT)

	
	def _delegateComputation(self, LinBox, sparse, parallel, minBlockSize, maxBlockSize, cores):
		if LinBox:
			low, high = self.complex.breaks[self.dimension], self.complex.breaks[self.dimension+1]
			
			def persist(filtration):
				essential = ComputePercolationEvents(
					self.matrices.full, filtration, self.dimension,
					self.complex.field, self.complex.breaks
				)

				essential = np.array(list(essential))
				essential.sort()
				
				return essential[(essential >= low) & (essential < high)]
			
		else:
			# If we can't/don't want to use LinBox, use inbuilt methods.
			Reducer = MatrixReduction(self.complex.field, parallel, minBlockSize, maxBlockSize, cores)
			Persistencer = Persistence(self.complex.field, self.complex.flattened, self.dimension)

			coboundary = np.zeros((self.cells, self.faces), dtype=FINT)
			rows = self.matrices.coboundary[::3]
			cols = self.matrices.coboundary[1::3]
			entries = (self.matrices.coboundary[2::3]%self.complex.field).astype(FINT)

			coboundary[rows,cols] = entries
			
			def persist(filtration):
				return Persistencer.TwistComputePercolationEvents(filtration)
			
		
		# If p == 2, then we want to use PHAT for persistence only. We have to some
		# unfortunately stupid computations though.
		if self.complex.field < 3:
			Persistencer = Persistence(self.complex.field, self.complex.flattened, self.dimension)
			t = self.complex.tranches

			# Indices of each cell and dimensions.
			dimensions = np.array(sum([
				[d]*(t[d,1]-t[d,0]) for d in range(len(t)) if d < self.dimension+2
			],[]), dtype=FINT)
			times = np.array(range(t[1][0], len(dimensions))).astype(FINT)

			def phattified(phatBoundary, dimensions, times, filtration):
				flattened = Persistencer.ReindexBoundary(filtration)
				
				dimensionalFlattened = [
					(d, sorted(f)) if d > 0 else (d, [])
					for d, f in zip(dimensions, flattened)
				]

				phatBoundary.columns = dimensionalFlattened
				_births, _deaths = zip(*phatBoundary.compute_persistence_pairs())
				births = set(_births)
				deaths = set(_deaths)

				return set(
					e for e in times-(births|deaths)
					if t[self.dimension][0] <= e < t[self.dimension][1]
				)

			def persist(filtration):
				return phattified(phat.boundary_matrix(), dimensions, set(times), filtration)

		self.persist = persist


	def filtrate(self):
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


	def initial(self):
		"""
		Computes an initial state for the model's Complex.

		Returns:
			A numpy `np.array` representing a vector of spin assignments.
		"""
		return self.RNG.integers(
			0, high=self.complex.field, dtype=FINT, size=self.faces
		)
	

	def proposal(self, time):
		"""
		Proposal scheme for generalized invaded-cluster evolution on the
		random-cluster model.

		Args:
			time (int): Step in the chain.

		Returns:
			A numpy array representing a vector of spin assignments.
		"""
		# Construct the filtration and find the essential cycles.
		filtration, shuffled = self.filtrate()
		essential = self.persist(filtration)

		j = 0
		low = self.complex.breaks[self.dimension]

		occupied = np.zeros((self.rank, self.nullity))

		for t in essential:
			occupiedIndices = shuffled[:t-low]
			occupied[j,occupiedIndices] = 1

		return occupied
	

	def assign(self, cocycle):
		"""
		Updates mappings from faces to spins and cubes to occupations.

		Args:
			cocycle (np.array): Cocycle on the subcomplex.
		
		Returns:
			Nothing.
		"""
		self.spins = cocycle

