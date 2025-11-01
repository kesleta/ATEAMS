
import numpy as np
import warnings

from ..arithmetic import ReducedKernelSample
from ..common import FINT, Matrices, TooSmallWarning, NumericalInstabilityWarning
from ..statistics import constant


class SwendsenWang():
	_name = "SwendsenWang"

	def __init__(
			self, C, dimension=1, field=2, temperature=constant(-0.6), initial=None,
			**kwargs
		):
		r"""
		Initializes Swendsen-Wang evolution on the Potts model.

		Args:
			C (Complex): The `Complex` object on which we'll be running experiments.
			dimension (int=1): The dimension of cells on which we're percolating.
			field (int=2): Field characteristic.
			temperature (Callable): A temperature schedule function which
				takes a single positive integer argument `t`, and returns the
				scheduled temperature at time `t`.
			initial (np.array): A vector of spin assignments to components.

		<center> <button type="button" class="collapsible" id="SwendsenWang-LanczosKernelSample-2"> Performance in \(\mathbb T^2_N\)</button> </center>
		.. include:: ./tables/SwendsenWang.LanczosKernelSample.2.html
		
		<center> <button type="button" class="collapsible" id="SwendsenWang-LanczosKernelSample-4"> Performance in \(\mathbb T^4_N\)</button> </center>
		.. include:: ./tables/SwendsenWang.LanczosKernelSample.4.html
		"""
		# Object access.
		self.complex = C
		self.temperature = temperature
		self.dimension = dimension
		self._returns = 2
		self.field = field

		# Force-recompute the matrices for a different dimension; creates
		# a set of orientations for fast elementwise products.
		self.matrices = Matrices()
		self.matrices.full = self.complex.matrices.full

		boundary, coboundary = self.complex.recomputeBoundaryMatrices(dimension)
		self.matrices.boundary = boundary
		self.matrices.coboundary = coboundary

		# Useful values to have later.
		self.cells = len(self.complex.Boundary[self.dimension])
		self.faces = len(self.complex.Boundary[self.dimension-1])

		# Check the dimensions of the boundary/coboundary matrices by comparing
		# the number of cells. LinBox is really sensitive to smaller-size matrices,
		# but can easily handle large ones.
		if self.cells*self.faces < 10000:
			warnings.warn(f"complex with {self.cells*self.faces} boundary matrix entries is too small for accurate matrix solves; may segfault.", TooSmallWarning, stacklevel=2)

		# Seed the random number generator.
		self.RNG = np.random.default_rng()

		# If no initial spin configuration is passed, initialize.
		if not initial: self.spins = self._initial()
		else: self.spins = (initial%self.field).astype(FINT)

		# Delegate computation.
		self._delegateComputation()

	
	def _delegateComputation(self):
		# If we use LinBox, keep everything as-is.
		def sample(zeros):
			# Not currently sure how to handle this... maybe we'll just "reject"
			# for now, come back, and sub something else in later. We shouldn't
			# be halting computation. For now, we should raise an exception that
			# the Chain catches, and warns the user by exiting with exit code
			# 1 or 2.
			try:
				return np.array(ReducedKernelSample(
					self.matrices.coboundary, zeros, 2*self.dimension,
					self.faces, self.field
				), dtype=FINT)
			except Exception as e:
				raise NumericalInstabilityWarning(e)
			
		self.sample = sample



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
		Proposal scheme for generalized Swendsen-Wang evolution on the Potts model.

		Args:
			time (int): Step in the chain.

		Returns:
			A numpy array representing a vector of spin assignments.
		"""
		# Compute the probability of choosing any individual cube in the complex.
		T = self.temperature(time)
		p = 1-np.exp(T)
		assert 0 <= p <= 1

		# Choose cubes to include.
		uniform = self.RNG.uniform(size=self.cells)
		include = np.nonzero(uniform < p)[0]
		boundary = self.complex.Boundary[self.dimension]
		q = self.field

		# Evaluate the current spin assignment (cochain).
		coefficients = self.spins[boundary]
		coefficients[:,1::2] = -coefficients[:,1::2]%q
		sums = coefficients.sum(axis=1)%q
		zeros = np.nonzero(sums == 0)[0]
		includedZeros = np.intersect1d(zeros, include)

		# Sample from the kernel of the coboundary matrix, and evaluate again
		spins = self.sample(includedZeros)

		# # Check whether the sample's correct.
		# coefficients = spins[boundary[includedZeros]]
		# coefficients[:,1::2] = -coefficients[:,1::2]%q
		# sums = coefficients.sum(axis=1)%q
		# print(np.nonzero(sums > 0)[0])

		# Create output vectors.
		occupied = np.zeros(self.cells, dtype=FINT)
		occupied[includedZeros] = 1

		return spins, occupied
	

	def _assign(self, cocycle):
		"""
		Updates mappings from faces to spins and cubes to occupations.

		Args:
			cocycle (np.array): Cocycle on the subcomplex.
		
		Returns:
			None.
		"""
		self.spins = cocycle

