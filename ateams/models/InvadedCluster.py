 
import numpy as np
from functools import partial
from math import comb

from ..arithmetic import (
	sampleFromKernel,
	boundaryMatrix,
	computePersistencePairs,
	evaluateCochain,
	reindexSparseBoundaryMatrix
)
from .Model import Model


class InvadedCluster(Model):
	name = "InvadedCluster"
	
	def __init__(self, L, homology=1, stop=lambda: 1):
		"""
		Initializes the plaquette invaded-cluster algorithm on the provided
		integer lattice, detecting percolation in the `homology`-th homology
		group.

		Args:
			L (Lattice): The `Lattice` object on which we'll be running experiments.
			homology (int=1): Computing the `homology`th homology group of the
				complex.
			stop (function): A function that returns the number of essential cycles
				found before sampling the next configuration.
		"""
		self.lattice = L
		self.homology = homology
		self.stop = stop

		# Set an initial spin configuration.
		self.spins = self.initial()

		# Specify the dimensions of each cell; time steps; non-randomized cells.
		dimensions = np.array(sum([
			[d]*(b-a) for d, (a,b) in self.lattice.tranches.items() if d < self.homology+2
		],[]))
		
		times = np.array(range(
			self.lattice.tranches[1][0], len(dimensions)
		))

		# Construct a "filtration template" of all the indices. At each step, we
		# determine which of the target indices correspond to satisfied `homology`
		# -dimensional cells, then shuffle the order of those indices.
		self.targetIndices = np.arange(*self.lattice.tranches[homology])
		self.zeroedTargetIndices = np.arange(len(self.lattice.boundary[homology]))

		# Set multiplicative inverses for the field we're working with.
		fieldInverses = dict(zip(
			[int(t) for t in self.lattice.field.Range(1, self.lattice.field.characteristic)],
			[int(t) for t in self.lattice.field.Range(1, self.lattice.field.characteristic)**(-1)]
		))

		# Create some pre-fab data structures to provide as fixed arguments to
		# the proposal method.
		premarked = list(range(self.lattice.tranches[0][1]))

		# Premake the "occupied cells" array; change the dimension of the lattice
		# to correspond to the provided dimension.
		self.rank = comb(len(self.lattice.corners), homology)
		self.nullity = len(self.lattice.boundary[homology])

		self.lattice.dimension = homology
		self.coboundary = boundaryMatrix(self.lattice.boundary, self.homology, self.lattice.field).T

		self.computeGiantCyclePairs = partial(
			computePersistencePairs,
			times,
			self.lattice.tranches[homology][1],
			premarked,
			dimensions,
			homology+1,
			self.lattice.tranches,
			{ d : len(v) for d, v in self.lattice.boundary.items() },
			fieldInverses,
			self.lattice.field.characteristic
		)
	
	
	def initial(self):
		"""
		Computes a uniform random initial spin configuration.

		Returns:
			A Galois `Array` of independent uniform random draws from the field
			of coefficients.
		"""
		return self.lattice.field.Random(len(self.lattice.boundary[self.homology-1]))
	

	def proposal(self, time):
		"""
		Proposal scheme for plaquette invaded-cluster. Each "step" in the chain
		is the set of spins from which a giant cycle first emerges (i.e. the
		"homological percolation" event).

		Args:
			time (int): Step in the chain; not used.

		Returns:
			A `(galois.FieldArray, numpy.ndarray, numpy.ndarray)` triplet representing
			the proposed spin configuration, the occupied plaquettes at each
			occurrence of homological percolation, and the satisfied plaquettes.
		"""
		boundary = self.lattice.boundary
		homology = self.homology
		cochain = self.spins

		# Determine the satisfied plaquettes and re-index the sparse boundary
		# matrix.
		cycles = evaluateCochain(boundary[homology], cochain)
		low, filtration, flattened, shuffledIndices, satisfiedIndices = reindexSparseBoundaryMatrix(
			cycles, boundary, self.lattice.flattened, homology, self.lattice.tranches,
			self.lattice.reindexed, self.targetIndices, self.zeroedTargetIndices
		)

		# Find essential cycles.
		essential = self.computeGiantCyclePairs(filtration, flattened)
	
		# Now, make sure we know when *all* the essential cycles are born.
		j = 0
		stop = self.stop()
		occupied = np.zeros((self.rank, self.nullity))
		satisfied = np.zeros(self.nullity)

		# For each essential cycle born...
		for t in sorted(essential):
			# Determine which cells were included at the time the cycle was born, and
			# construct three assignments: the *occupied* cells, the *satisfied* cells,
			# and a new spin assignment on the *faces* of the cells.
			occupiedIndices = shuffledIndices[:t-low]
			occupied[j][occupiedIndices] = 1

			# Only sample the next cocycle from the time we homologically percolate,
			# not after.
			if (j+1) == stop:
				spins = sampleFromKernel(self.coboundary, self.lattice.field, includes=occupiedIndices)

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
