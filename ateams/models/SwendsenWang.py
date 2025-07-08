
import numpy as np
import warnings

from ..arithmetic import LanczosKernelSample, KernelSample, MatrixReduction
from ..common import MINT, FINT, Matrices, TooSmallWarning, NumericalInstabilityWarning
from ..statistics import constant
from .Model import Model


class SwendsenWang():
	_name = "SwendsenWang"

	def __init__(
			self, C, dimension=1, field=2, temperature=constant(-0.6), initial=None, LinBox=True,
			parallel=False, minBlockSize=32, maxBlockSize=64, cores=4,
			maxTries=16
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
			LinBox (bool=True): Uses fast LinBox routines instead of slow inbuilt
				ones. WARNING: using inbuilt methods may dramatically increase
				computation time.
			parallel (bool=False): Should matrix computations be done in parallel? (Uses C/C++).
			minBlockSize (int=32): If `parallel` is truthy, this is the smallest
				number of columns processed in parallel.
			maxBlockSize (int=64): If `parallel` is truthy, this is the largest
				number of columns processed in parallel.
			cores (int=4): Number of available CPUs/cores/threads on the machine.
			maxTries (int=16): The number of attempts LinBox makes to sample a nonzero
				vector in the kernel of the coboundary matrix.


		If the first vector sampled
		by the Lanczos algorithm is all zeros, we perform the following
		steps until a nonzero vector is found, or we exhaust the number
		of attempts:

		1. sample once from \(A\) with no preconditioner;
		2. sample up to twice from \(D_0A\), where \(D_0\) is a random diagonal matrix;
		3. sample up to twice from \(A^\top D_1 A\), where \(D_1\) is a random diagonal matrix;
		4. sample for the remainder of the attempts from \(D_2 A^\top D_3 A D_2\), where \(D_2, D_3\) are random diagonal matrices.
			
		If we spend the entire budget of attempts, it is likely that the
		result returned is either the all-zeros vector.
		
		Included below are performance statistics for various configurations of
		`SwendsenWang`. Each configuration completed 100 iterations on Pangolin,
		a Dell Precision 5280 workstation with an 18-core@1.305GhZ Intel Xeon W-2295
		CPU. **Warns** indicates the number of `NumericalInstabilityWarning`s caught
		during the run, and **Zeros** the number of all-zeros vectors returned.

		</br>
		</br>
		<center> <strong> Samples in \(\mathbb T^2_N\) </strong> </center>
		.. include:: ./tables/SwendsenWang.LanczosKernelSample.2.html

		</br>
		</br>
		<center> <strong> Samples in \(\mathbb T^4_N\) </strong> </center>
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
		if self.cells*self.faces < 10000 and LinBox:
			warnings.warn(f"complex with {self.cells*self.faces} boundary matrix entries is too small for accurate matrix solves; may segfault.", TooSmallWarning, stacklevel=2)

		# Seed the random number generator.
		self.RNG = np.random.default_rng()

		# If no initial spin configuration is passed, initialize.
		if not initial: self.spins = self._initial()
		else: self.spins = (initial%self.field).astype(FINT)

		# Delegate computation.
		self._delegateComputation(LinBox, parallel, minBlockSize, maxBlockSize, cores, maxTries)

	
	def _delegateComputation(self, LinBox, parallel, minBlockSize, maxBlockSize, cores, maxTries):
		# If we use LinBox, keep everything as-is.
		if LinBox:
			def sample(zeros):
				# Not currently sure how to handle this... maybe we'll just "reject"
				# for now, come back, and sub something else in later. We shouldn't
				# be halting computation. For now, we should raise an exception that
				# the Chain catches, and warns the user by exiting with exit code
				# 1 or 2.
				try:
					return np.array(LanczosKernelSample(
						self.matrices.coboundary, zeros, 2*self.dimension,
						self.faces, self.field, maxTries=maxTries
					), dtype=FINT)
				except Exception as e:
					raise NumericalInstabilityWarning(e)

		
		# If we use inbuilt routines, we need to actually construct the matrix
		# form of the coboundary matrix, which is extremely sparse (and really big).
		# Makes things slow.
		else:
			coboundary = np.zeros((self.cells, self.faces), dtype=FINT)
			rows = self.matrices.coboundary[::3]
			cols = self.matrices.coboundary[1::3]
			entries = (self.matrices.coboundary[2::3]%self.field).astype(FINT)

			coboundary[rows,cols] = entries
			Reducer = MatrixReduction(self.field, parallel, minBlockSize, maxBlockSize, cores)

			def sample(zeros):
				return KernelSample(Reducer, coboundary.take(zeros, axis=0)).astype(FINT)
			
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
		coefficients = self.spins[boundary[include]]
		coefficients[:,1::2] = -coefficients[:,1::2]%q
		sums = coefficients.sum(axis=1)%q
		zeros = np.nonzero(sums == 0)[0]

		# Create output vectors.
		satisfied = np.zeros(self.cells, dtype=FINT)
		satisfied[zeros] = 1

		# Sample from the kernel of the coboundary matrix.
		spins = self.sample(zeros)

		return spins, satisfied
	

	def _assign(self, cocycle):
		"""
		Updates mappings from faces to spins and cubes to occupations.

		Args:
			cocycle (np.array): Cocycle on the subcomplex.
		
		Returns:
			None.
		"""
		self.spins = cocycle

