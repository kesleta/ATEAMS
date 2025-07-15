
import numpy as np
import warnings

from math import comb

from ..arithmetic import (
	Twist, LanczosKernelSample, MatrixReduction,
	Persistence, KernelSample, ComputePersistencePairs
)
from ..common import Matrices, TooSmallWarning, NumericalInstabilityWarning, FINT


class InvadedCluster():
	_name = "InvadedCluster"
	
	def __init__(
			self, C, dimension=1, field=2, initial=None, stop=lambda: 1, LinBox=True, maxTries=16,
			parallel=False, minBlockSize=32, maxBlockSize=64, cores=4, **kwargs
		):
		"""
		Initializes the plaquette invaded-cluster algorithm on the provided
		integer complex, detecting percolation in the `homology`-th homology
		group.

		Args:
			C (Complex): The `Complex` object on which we'll be running experiments.
			dimension (int=1): The dimension of cells on which we're percolating.
			field (int=2): Field characteristic.
			initial (np.array): A vector of spin assignments to components.
			stop (function): A function that returns the number of essential cycles
				found before sampling the next configuration.
			LinBox (bool=True): Uses fast LinBox routines instead of slow inbuilt
				ones. WARNING: using inbuilt methods may dramatically increase
				computation time.
			maxTries (int=16): The number of attempts LinBox makes to sample a nonzero
				vector in the kernel of the coboundary matrix \(A\). See more
				discussion in the `SwendsenWang` model.
			parallel (boolean): Should matrix computations be done in parallel? (Uses C/C++).
			minBlockSize (int=32): If `parallel` is truthy, this is the smallest
				number of columns processed in parallel.
			maxBlockSize (int=64): If `parallel` is truthy, this is the largest
				number of columns processed in parallel.
			cores (int=4): Number of available CPUs/cores/threads on the machine.

		<center> <button type="button" class="collapsible" id="InvadedCluster-Persistence-2">Performance in \(\mathbb T^2_N\)</button> </center>
		..include:: ./tables/InvadedCluster.Persistence.2.html

		<center> <button type="button" class="collapsible" id="InvadedCluster-Persistence-4">Performance in \(\mathbb T^4_N\)</button> </center>
		..include:: ./tables/InvadedCluster.Persistence.4.html
		"""
		# Object access.
		self.complex = C
		self.dimension = dimension
		self.stop = stop
		self._returns = 3
		self.field = field


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


		# Check the dimensions of the boundary/coboundary matrices by comparing
		# the number of cells. LinBox is really sensitive to smaller-size matrices,
		# but can easily handle large ones.
		if self.cells*self.faces < 10000 and LinBox:
			warnings.warn(f"complex with {self.cells*self.faces} boundary matrix entries is too small for accurate matrix solves; may segfault.", TooSmallWarning, stacklevel=2)


		# Premake the "occupied cells" array; change the dimension of the complex
		# to correspond to the provided dimension.
		self.rank = comb(len(self.complex.corners), self.dimension)
		self.nullity = len(self.complex.Boundary[self.dimension])


		# Delegates computation for persistence and cocycle sampling.
		self._delegateComputation(LinBox, parallel, minBlockSize, maxBlockSize, cores, maxTries)


		# Seed the random number generator.
		self.RNG = np.random.default_rng()


		# If no initial spin configuration is passed, initialize.
		if not initial: self.spins = self._initial()
		else: self.spins = (initial%self.field).astype(FINT)

	
	def _delegateComputation(self, LinBox, parallel, minBlockSize, maxBlockSize, cores, maxTries):
		low, high = self.complex.breaks[self.dimension], self.complex.breaks[self.dimension+1]

		# If we're using LinBox, our sampling method is Lanczos regardless of
		# dimension.
		if LinBox:
			def sample(zeros):
				try:
					return np.array(LanczosKernelSample(
						self.matrices.coboundary, zeros, 2*self.dimension,
						self.faces, self.field, maxTries
					), dtype=FINT)
				except Exception as e:
					raise NumericalInstabilityWarning(e)
		
		# If we're using LinBox and the characteristic of our field is greater
		# than 2, we use the twist_reduce variant implemented in this library.
		if LinBox and self.field > 2:
			Twister = Twist(self.field, self.matrices.full, self.complex.breaks, self.cellCount, self.dimension)

			def persist(filtration):
				essential = Twister.LinearComputePercolationEvents(filtration)
				# essential = Twister.ComputePercolationEvents(filtration)
				# essential = Twister.ZpComputePercolationEvents(filtration)
				essential = np.array(list(essential))
				essential.sort()
				
				return essential[(essential >= low) & (essential < high)]
		
		# If we're using LinBox and the field we're computing over *is* two,
		# use PHAT.
		elif LinBox and self.field < 3:
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
				essential = ComputePersistencePairs(self.matrices.full, filtration, self.dimension, self.complex.breaks)
				return whittle(essential)
			
		if not LinBox:
			# If we can't/don't want to use LinBox, use inbuilt methods.
			Reducer = MatrixReduction(self.field, parallel, minBlockSize, maxBlockSize, cores)
			Persistencer = Persistence(self.field, self.complex.flattened, self.dimension)

			coboundary = np.zeros((self.cells, self.faces), dtype=FINT)
			rows = self.matrices.coboundary[::3]
			cols = self.matrices.coboundary[1::3]
			entries = (self.matrices.coboundary[2::3]%self.field).astype(FINT)

			coboundary[rows,cols] = entries
			
			def persist(filtration):
				return Persistencer.TwistComputePercolationEvents(filtration)
			
			def sample(zeros):
				return KernelSample(Reducer, coboundary.take(zeros, axis=0)).astype(FINT)
		

		self.sample = sample
		self.persist = persist


	def _filtrate(self, cochain):
		"""
		Constructs a filtration based on the evaluation of the cochain.
		"""
		# Find which cubes get zeroed out (i.e. are sent to zero by the cocycle).
		boundary = self.complex.Boundary[self.dimension]
		q = self.field
		
		coefficients = cochain[boundary]
		coefficients[:,1::2] = -coefficients[:,1::2]%q
		sums = coefficients.sum(axis=1)%q

		# POSSIBLY INEFFICIENT!! TAKE A LOOK DUMMY
		satisfied = np.nonzero(sums==0)[0]
		unsatisfied = np.nonzero(sums>0)[0]
		m = satisfied.shape[0]

		# Construct the filtration.
		filtration = np.arange(self.cellCount)
		low = self.complex.breaks[self.dimension]
		high = self.complex.breaks[self.dimension+1]

		shuffled = np.random.permutation(satisfied)
		filtration[low:low+m] = self.target[shuffled]
		filtration[low+m:high] = self.target[unsatisfied]

		# with open("filtration.txt", "w") as w:
		# 	for f in np.arange(self.cellCount): w.write(f"{f}\n")

		# with open("bd.txt", "w") as w:
		# 	for f in self.matrices.full: w.write(f"{f}\n")

		return filtration, shuffled, satisfied


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
		filtration, shuffledIndices, satisfiedIndices = self._filtrate(self.spins)
		essential = self.persist(filtration)

		j = 0
		low = self.complex.breaks[self.dimension]

		stop = self.stop()
		occupied = np.zeros((self.rank, self.nullity))
		satisfied = np.zeros(self.nullity)

		for t in essential:
			occupiedIndices = shuffledIndices[:t-low]
			occupied[j,occupiedIndices] = 1

			if (j+1) == stop: spins = self.sample(occupiedIndices)

			j += 1

		satisfied[satisfiedIndices] = 1

		return spins, occupied, satisfied
	

	def _assign(self, cocycle):
		"""
		Updates mappings from faces to spins and cubes to occupations.

		Args:
			cocycle (np.array): Cocycle on the subcomplex.
		
		Returns:
			Nothing.
		"""
		self.spins = cocycle

