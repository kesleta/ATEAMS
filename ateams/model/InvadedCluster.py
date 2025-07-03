
import numpy as np
from math import comb

from ..arithmetic import (
	ComputePercolationEvents, LanczosKernelSample, MatrixReduction,
	Persistence, FINT
)
from ..common import Matrices
from .Model import Model


class InvadedCluster(Model):
	name = "InvadedCluster"
	
	def __init__(
			self, C, homology=1, initial=None, stop=lambda: 1, LinBox=True, sparse=True,
			parallel=False, minBlockSize=32, maxBlockSize=64, cores=4
		):
		"""
		Initializes the plaquette invaded-cluster algorithm on the provided
		integer lattice, detecting percolation in the `homology`-th homology
		group.

		Args:
			C (Complex): The `Complex` object on which we'll be running experiments.
			homology (int=1): Computing the `homology`th homology group of the
				complex.
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
		self.homology = homology
		self.stop = stop

		# Force-recompute the matrices for a different dimension; creates
		# a set of orientations for fast elementwise products.
		self.matrices = Matrices()
		self.matrices.full = self.complex.matrices.full

		boundary, coboundary = self.complex.recomputeBoundaryMatrices(homology)
		self.matrices.boundary = boundary
		self.matrices.coboundary = coboundary


		# Useful values to have later.
		self.cellCount = len(self.complex.flattened)
		self.cells = len(self.complex.Boundary[self.homology])
		self.faces = len(self.complex.Boundary[self.homology-1])
		self.orientations = np.tile([-1,1], self.homology).astype(FINT)

		self.target = np.arange(self.complex.breaks[self.homology], self.complex.breaks[self.homology+1])

		# Premake the "occupied cells" array; change the dimension of the lattice
		# to correspond to the provided dimension.
		self.rank = comb(len(self.complex.corners), self.homology)
		self.nullity = len(self.complex.Boundary[self.homology])

		# Delegates computation for persistence and cocycle sampling.
		self._delegateComputation(LinBox, sparse, parallel, minBlockSize, maxBlockSize, cores)

		# Seed the random number generator.
		self.RNG = np.random.default_rng()

		# If no initial spin configuration is passed, initialize.
		if not initial: self.spins = self.initial()
		else: self.spins = (initial%self.complex.field).astype(FINT)

	
	def _delegateComputation(self, LinBox, sparse, parallel, minBlockSize, maxBlockSize, cores):
		if LinBox:
			low, high = self.complex.breaks[self.homology], self.complex.breaks[self.homology+1]

			def sample(zeros):
				return np.array(LanczosKernelSample(
					self.matrices.coboundary, zeros, 2*self.homology,
					self.faces, self.complex.field
				), dtype=FINT)
			
			def persist(filtration):
				# Get the "essential" cycles, then whittle down the ones we care
				# about. We want to return the essential cycles in sorted order.
				essential = ComputePercolationEvents(
					self.complex.matrices.full, filtration, self.homology,
					self.complex.field, self.complex.breaks
				)

				essential = np.array(list(essential))
				essential.sort()
				
				return essential[(essential >= low) & (essential < high)]

		self.sample = sample
		self.persist = persist


	def filtrate(self, cochain):
		"""
		Constructs a filtration based on the evaluation of the cochain.
		"""
		# Find which cubes get zeroed out (i.e. are sent to zero by the cocycle).
		boundary = self.complex.Boundary[self.homology]
		q = self.complex.field

		coefficients = (cochain[boundary]*self.orientations)%q
		sums = coefficients.sum(axis=1)%q

		# POSSIBLY INEFFICIENT!! TAKE A LOOK DUMMY
		satisfied = np.nonzero(sums==0)[0]
		unsatisfied = np.nonzero(sums>0)[0]
		m = satisfied.shape[0]

		# Construct the filtration.
		filtration = np.arange(self.cellCount)
		low = self.complex.breaks[self.homology]
		high = self.complex.breaks[self.homology+1]

		shuffled = np.random.permutation(satisfied)
		filtration[low:low+m] = self.target[shuffled]
		filtration[low+m:high] = self.target[unsatisfied]

		return filtration, shuffled, satisfied
		return np.arange(self.cellCount), shuffled, satisfied


	def initial(self):
		"""
		Computes an initial state for the model's Lattice.

		Returns:
			A Galois `Array` representing a vector of spin assignments.
		"""
		return self.RNG.integers(
			0, high=self.complex.field, dtype=FINT, size=self.faces
		)
	

	def proposal(self, time):
		# Construct the filtration and find the essential cycles.
		filtration, shuffledIndices, satisfiedIndices = self.filtrate(self.spins)
		essential = self.persist(filtration)

		j = 0
		low = self.complex.breaks[self.homology]

		stop = self.stop()
		occupied = np.zeros((self.rank, self.nullity))
		satisfied = np.zeros(self.nullity)

		for t in essential:
			occupiedIndices = shuffledIndices[:t-low]
			occupied[j,occupiedIndices] = 1

			if (j+1) == stop:
				spins = self.sample(occupiedIndices)

			j += 1

		satisfied[satisfiedIndices] = 1

		return spins, occupied, satisfied
	

	def assign(self, cocycle):
		"""
		Updates mappings from faces to spins and cubes to occupations.

		Args:
			cocycle (galois.FieldArray): Cocycle on the sublattice.
		
		Returns:
			Nothing.
		"""
		self.spins = cocycle

