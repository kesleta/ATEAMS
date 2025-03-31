 
import numpy as np
import phat
from functools import partial
from math import comb

from ..arithmetic import (
	SparseSampleFromKernel,
	SampleFromKernel,
	sampleFromKernel,
	boundaryMatrix,
	KernelBasis,
	computeGiantCyclePairs,
	evaluateCochain,
	reindexSparseBoundaryMatrix
)
from .Model import Model

DTYPE = np.int64


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

		# Set an initial spin configuration. After we assign this, we use PHAT
		# if we're working over Z/2Z coefficients, or the Zomorodian/Edelsbrunner
		# algorithm otherwise.
		self.spins = self.initial()
		self._initializeProposal()
		self.coboundary = boundaryMatrix(L.boundary, self.homology, L.field).T

		# Pre-load arguments for re-indexing the boundary matrix.
		self.reindexSparseBoundaryMatrix = partial(
			reindexSparseBoundaryMatrix,
			self.lattice.boundary[homology+1].astype(DTYPE),
			self.lattice.flattened,
			int(self.homology),
			self.lattice.tranches.astype(DTYPE),
			self.lattice.reindexed[homology].astype(DTYPE),
			self.targetIndices.astype(DTYPE),
			self.zeroedTargetIndices.astype(DTYPE)
		)

	
	def _initializeProposal(self):
		t = self.lattice.tranches
		p = self.lattice.field.characteristic

		# Premake the "occupied cells" array; change the dimension of the lattice
		# to correspond to the provided dimension.
		self.rank = comb(len(self.lattice.corners), self.homology)
		self.nullity = len(self.lattice.boundary[self.homology])

		# Construct a "filtration template" of all the indices. At each step, we
		# determine which of the target indices correspond to satisfied `homology`
		# -dimensional cells, then shuffle the order of those indices.
		self.targetIndices = np.arange(*t[self.homology], dtype=DTYPE)
		self.zeroedTargetIndices = np.arange(len(self.lattice.boundary[self.homology]), dtype=DTYPE)

		# Indices of each cell and dimensions.
		dimensions = np.array(sum([
			[d]*(t[d,1]-t[d,0]) for d in range(len(t)) if d < self.homology+2
		],[]), dtype=DTYPE)

		times = np.array(range(t[1][0], len(dimensions))).astype(DTYPE)

		# If p > 2, we use the simplified Edelsbrunner/Zomorodian persistence
		# algorithm; otherwise, we use PHAT (Bauer et al.).
		if p > 2:
			# Find the max over the dimensions and specify a "blank" array the max
			# width.
			zeros = np.zeros((2, t[:,1].max()//8), dtype=DTYPE)

			# Set multiplicative inverses for the field we're working with.
			fieldInverses = np.array([0] + list(self.lattice.field.Range(1, p)**(-1))).astype(DTYPE)

			# Create some pre-fab data structures to provide as fixed arguments to
			# the proposal method.
			low = t[0][1]
			premarked = np.array(list(range(low)), dtype=DTYPE)

			self.lattice.dimension = self.homology
			self.coboundary = boundaryMatrix(self.lattice.boundary, self.homology, self.lattice.field).T

			# Create arithmetic lookup tables.
			addition = np.zeros((p,p), dtype=DTYPE)
			for j in range(p): addition[:,j] = (np.arange(p, dtype=DTYPE)+j)%p

			subtraction = np.zeros((p,p), dtype=DTYPE)
			for j in range(p): subtraction[:,j] = (np.arange(p, dtype=DTYPE)-j)%p

			multiplication = np.zeros((p,p), dtype=DTYPE)
			for j in range(p): multiplication[:,j] = (np.arange(p, dtype=DTYPE)*j)%p

			powers = np.full(zeros.shape[1], -1, dtype=DTYPE)
			powers[1::2] = -powers[1::2]
			powers = powers%p

			self.computeGiantCyclePairs = partial(
				computeGiantCyclePairs,
				times,
				premarked,
				dimensions,
				t,
				fieldInverses,
				p,
				self.homology+1,
				t[self.homology+1][1],
				np.empty((2,0), dtype=DTYPE),
				zeros,
				np.empty(zeros.shape[1], dtype=DTYPE),
				addition.astype(DTYPE),
				subtraction.astype(DTYPE),
				multiplication.astype(DTYPE),
				powers.astype(DTYPE)
			)
		else:
			def phattified(phatBoundary, dimensions, times, filtration, flattened, zeros):
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
					if t[self.homology][0] <= e < t[self.homology][1]
				)

			self.computeGiantCyclePairs = partial(
				phattified,
				phat.boundary_matrix(),
				dimensions,
				set(times)
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
		cycles = evaluateCochain(boundary[homology], cochain).astype(DTYPE)
		low, filtration, flattened, shuffledIndices, satisfiedIndices = self.reindexSparseBoundaryMatrix(cycles)

		# Find essential cycles.
		essential = self.computeGiantCyclePairs(
			filtration,
			flattened,
			np.zeros(filtration.shape[0], dtype=DTYPE)
		)
	
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


class CInvadedCluster(Model):
	name = "InvadedCluster"
	
	def __init__(self, L, homology=1, stop=lambda: 1, parallel=False):
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
			parallel (boolean): Should matrix computations be done in parallel?
		"""
		self.lattice = L
		self.homology = homology
		self.stop = stop
		self.parallel = parallel

		# Set an initial spin configuration. After we assign this, we use PHAT
		# if we're working over Z/2Z coefficients, or the Zomorodian/Edelsbrunner
		# algorithm otherwise.
		self.spins = self.initial()
		self._initializeProposal()
		self.coboundary = boundaryMatrix(L.boundary, self.homology, L.field).T

		# Pre-load arguments for re-indexing the boundary matrix.
		self.reindexSparseBoundaryMatrix = partial(
			reindexSparseBoundaryMatrix,
			self.lattice.boundary[homology+1].astype(DTYPE),
			self.lattice.flattened,
			int(self.homology),
			self.lattice.tranches.astype(DTYPE),
			self.lattice.reindexed[homology].astype(DTYPE),
			self.targetIndices.astype(DTYPE),
			self.zeroedTargetIndices.astype(DTYPE)
		)

	
	def _initializeProposal(self):
		t = self.lattice.tranches
		p = self.lattice.field.characteristic

		# Premake the "occupied cells" array; change the dimension of the lattice
		# to correspond to the provided dimension.
		self.rank = comb(len(self.lattice.corners), self.homology)
		self.nullity = len(self.lattice.boundary[self.homology])

		# Construct a "filtration template" of all the indices. At each step, we
		# determine which of the target indices correspond to satisfied `homology`
		# -dimensional cells, then shuffle the order of those indices.
		self.targetIndices = np.arange(*t[self.homology], dtype=DTYPE)
		self.zeroedTargetIndices = np.arange(len(self.lattice.boundary[self.homology]), dtype=DTYPE)

		# Indices of each cell and dimensions.
		dimensions = np.array(sum([
			[d]*(t[d,1]-t[d,0]) for d in range(len(t)) if d < self.homology+2
		],[]), dtype=DTYPE)

		times = np.array(range(t[1][0], len(dimensions))).astype(DTYPE)

		# Find the max over the dimensions and specify a "blank" array the max
		# width.
		zeros = np.zeros((2, t[:,1].max()//8), dtype=DTYPE)

		# Set multiplicative inverses for the field we're working with.
		fieldInverses = np.array([0] + list(self.lattice.field.Range(1, p)**(-1))).astype(DTYPE)

		# Create some pre-fab data structures to provide as fixed arguments to
		# the proposal method.
		low = t[0][1]
		premarked = np.array(list(range(low)), dtype=DTYPE)

		self.lattice.dimension = self.homology
		self.coboundary = boundaryMatrix(self.lattice.boundary, self.homology, self.lattice.field).T
		self.identity = np.identity(self.coboundary.shape[1], dtype=DTYPE)

		# Create arithmetic lookup tables.
		addition = np.zeros((p,p), dtype=DTYPE)
		for j in range(p): addition[:,j] = (np.arange(p, dtype=DTYPE)+j)%p

		subtraction = np.zeros((p,p), dtype=DTYPE)
		for j in range(p): subtraction[:,j] = (np.arange(p, dtype=DTYPE)-j)%p

		negation = np.array([-q%p for q in range(0, p)], dtype=DTYPE)
		inverses = np.array([0] + list(self.lattice.field.Range(1, p)**(-1))).astype(DTYPE)

		multiplication = np.zeros((p,p), dtype=DTYPE)
		for j in range(p): multiplication[:,j] = (np.arange(p, dtype=DTYPE)*j)%p

		powers = np.full(zeros.shape[1], -1, dtype=DTYPE)
		powers[1::2] = -powers[1::2]
		powers = powers%p

		pivots = np.zeros(self.coboundary.shape[1], dtype=DTYPE)
		result = np.zeros(self.coboundary.shape[1], dtype=DTYPE)
		store = np.zeros(self.coboundary.shape[1], dtype=DTYPE)

		# If p > 2, we use the simplified Edelsbrunner/Zomorodian persistence
		# algorithm; otherwise, we use PHAT (Bauer et al.).
		if p > 2:
			self.computeGiantCyclePairs = partial(
				computeGiantCyclePairs,
				times,
				premarked,
				dimensions,
				t,
				fieldInverses,
				p,
				self.homology+1,
				t[self.homology+1][1],
				np.empty((2,0), dtype=DTYPE),
				zeros,
				np.empty(zeros.shape[1], dtype=DTYPE),
				addition.astype(DTYPE),
				subtraction.astype(DTYPE),
				multiplication.astype(DTYPE),
				powers.astype(DTYPE)
			)
		else:
			def phattified(phatBoundary, dimensions, times, filtration, flattened, zeros):
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
					if t[self.homology][0] <= e < t[self.homology][1]
				)

			self.computeGiantCyclePairs = partial(
				phattified,
				phat.boundary_matrix(),
				dimensions,
				set(times)
			)

		self.SampleFromKernel = partial(
			SparseSampleFromKernel,
			p,
			pivots,
			store,
			result,
			addition,
			subtraction,
			negation,
			multiplication,
			inverses,
			self.parallel
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
		cycles = evaluateCochain(boundary[homology], cochain).astype(DTYPE)
		low, filtration, flattened, shuffledIndices, satisfiedIndices = self.reindexSparseBoundaryMatrix(cycles)

		# Find essential cycles.
		essential = self.computeGiantCyclePairs(
			filtration,
			flattened,
			np.zeros(filtration.shape[0], dtype=DTYPE)
		)
	
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
				subbasis = np.concatenate([
					self.coboundary.take(occupiedIndices, axis=0).T,
					self.identity
				], axis=1, dtype=DTYPE)
				spins = self.SampleFromKernel(subbasis, len(occupiedIndices))
				spins = self.lattice.field(spins)

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