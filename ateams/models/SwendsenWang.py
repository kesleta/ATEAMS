
import numpy as np
import warnings

from ..arithmetic import LanczosKernelSample, KernelSample, MatrixReduction
from ..common import MINT, FINT, Matrices, TooSmallWarning
from ..statistics import constant
from .Model import Model


class SwendsenWang(Model):
	name = "SwendsenWang"

	def __init__(
			self, C, dimension=1, temperature=constant(-0.6), initial=None, LinBox=True,
			sparse=True, parallel=False, minBlockSize=32, maxBlockSize=64, cores=4
		):
		"""
		Initializes Swendsen-Wang evolution on the Potts model.

		Args:
			C (Complex): The `Complex` object on which we'll be running experiments.
			dimension (int=1): The dimension of cells on which we're percolating.
			temperature (Callable): A temperature schedule function which
				takes a single positive integer argument `t`, and returns the
				scheduled temperature at time `t`.
			initial (np.array): A vector of spin assignments to components.
			LinBox (bool=True): Uses fast LinBox routines instead of slow inbuilt
				ones. WARNING: using inbuilt methods may dramatically increase
				computation time.
			sparse (bool=True): Should matrices be formatted sparsely? (Uses C/C++).
			parallel (bool=False): Should matrix computations be done in parallel? (Uses C/C++).
			minBlockSize (int=32): If `parallel` is truthy, this is the smallest
				number of columns processed in parallel.
			maxBlockSize (int=64): If `parallel` is truthy, this is the largest
				number of columns processed in parallel.
			cores (int=4): Number of available CPUs/cores/threads on the machine.
		"""
		# Object access.
		self.complex = C
		self.temperature = temperature
		self.dimension = dimension


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
		self.orientations = np.tile([-1,1], self.dimension).astype(FINT)


		# Check the dimensions of the boundary/coboundary matrices by comparing
		# the number of cells. LinBox is really sensitive to smaller-size matrices,
		# but can easily handle large ones.
		if self.cells*self.faces < 10000 and LinBox:
			warnings.warn(f"complex with {self.cells*self.faces} boundary matrix entries is too small for accurate matrix solves; may segfault.", TooSmallWarning, stacklevel=2)


		# Seed the random number generator.
		self.RNG = np.random.default_rng()

		# If no initial spin configuration is passed, initialize.
		if not initial: self.spins = self.initial()
		else: self.spins = (initial%self.complex.field).astype(FINT)

		# Delegate computation.
		self._delegateComputation(LinBox, sparse, parallel, minBlockSize, maxBlockSize, cores)

	
	def _delegateComputation(self, LinBox, sparse, parallel, minBlockSize, maxBlockSize, cores):
		# If we use LinBox, keep everything as-is.
		if LinBox:
			def sample(zeros):
				if zeros.shape[0] < 1: return self.spins;
				
				return np.array(LanczosKernelSample(
					self.matrices.coboundary, zeros, 2*self.dimension,
					self.faces, self.complex.field
				), dtype=FINT)
		
		# If we use inbuilt routines, we need to actually construct the matrix
		# form of the coboundary matrix, which is extremely sparse (and really big).
		# Makes things slow.
		else:
			coboundary = np.zeros((self.cells, self.faces), dtype=FINT)
			rows = self.matrices.coboundary[::3]
			cols = self.matrices.coboundary[1::3]
			entries = (self.matrices.coboundary[2::3]%self.complex.field).astype(FINT)

			coboundary[rows,cols] = entries
			Reducer = MatrixReduction(self.complex.field, parallel, minBlockSize, maxBlockSize, cores)

			def sample(zeros):
				return KernelSample(Reducer, coboundary.take(zeros, axis=0)).astype(FINT)
			
		self.sample = sample



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
		q = self.complex.field

		# Evaluate the current spin assignment (cochain).
		coefficients = (self.spins[boundary[include]]*self.orientations)%q
		sums = coefficients.sum(axis=1)%q
		zeros = np.nonzero(sums == 0)[0]

		# Create output vectors.
		satisfied = np.zeros(self.cells, dtype=FINT)
		satisfied[zeros] = 1

		# Sample from the kernel of the coboundary matrix.
		spins = self.sample(zeros)

		return spins, satisfied
	

	def assign(self, cocycle):
		"""
		Updates mappings from faces to spins and cubes to occupations.

		Args:
			cocycle (galois.FieldArray): Cocycle on the sublattice.
		
		Returns:
			None.
		"""
		self.spins = cocycle

