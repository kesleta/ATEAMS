
import numpy as np
import warnings

from ..arithmetic import SubLanczosKernelSample
from ..common import FINT, TooSmallWarning, Matrices, Bunch, NumericalInstabilityWarning


class Nienhuis():
	_name = "Nienhuis"

	def __init__(
			self, C, q1, q2, dimension=2, field=2, initial=None, maxTries=16,
			**kwargs
		):
		"""
		Initializes the self-dual Nienhuis model.

		Args:
			C (Complex): The `Complex` object on which we'll be running experiments.
				Here, the dimension is assumed to be 2.
			q1 (float): Edge coupling parameter.
			q2 (float): Face (i.e. cell) coupling parameter.
			dimension (int): Dimension targeted by the model.
			field (int=2): Field characteristic.
			initial (np.array): A vector of spin assignments to components.
			maxTries (int=16): The number of attempts LinBox makes to sample a nonzero
				vector in the kernel of the coboundary matrix.
		
		<center> <button type="button" class="collapsible" id="Nienhuis-SubLanczosKernelSample-2"> Performance in \(\mathbb T^2_N\)</button> </center>
		..include:: ./tables/Nienhuis.SubLanczosKernelSample.2.html
		"""
		# Object access.
		self.complex = C
		self.dimension = dimension
		self._returns = 3
		self.field = field

		# Force-recompute the matrices for a different dimension; creates
		# a set of orientations for fast elementwise products. Here, the dimension
		# is implicitly assumed to be 2.
		self.matrices = Matrices()
		self.matrices.full = self.complex.matrices.full

		boundary, coboundary = self.complex.recomputeBoundaryMatrices(self.dimension)
		self.matrices.boundary = boundary
		self.matrices.coboundary = coboundary

		# Useful values to have later.
		self.cells = len(self.complex.Boundary[self.dimension])
		self.faces = len(self.complex.Boundary[self.dimension-1])
		self.faceIndices = np.arange(self.faces)

		self.orientations = Bunch()
		self.orientations.cells = np.tile([-1,1], self.dimension).astype(FINT)
		self.orientations.faces = np.tile([-1,1], self.dimension-1).astype(FINT)

		# Coupling parameters.
		self.coupling = Bunch()
		self.coupling.face = q1
		self.coupling.cell = q2

		self.prob = Bunch()
		self.prob.face = self.coupling.face/(1+self.coupling.face)
		self.prob.cell = self.coupling.cell/(1+self.coupling.cell)


		# Check the dimensions of the boundary/coboundary matrices by comparing
		# the number of cells. LinBox is really sensitive to smaller-size matrices,
		# but can easily handle large ones.
		if self.cells*self.faces < 100000:
			warnings.warn(f"complex with {self.cells*self.faces} boundary matrix entries is too small for accurate matrix solves; may segfault.", TooSmallWarning, stacklevel=2)

		# Seed the random number generator.
		self.RNG = np.random.default_rng()

		# If no initial spin configuration is passed, initialize.
		if not initial: self.spins = self._initial()
		else: self.spins = (initial%self.field).astype(FINT)

		# Delegate computation.
		self._delegateComputation(maxTries)


	def _delegateComputation(self, maxTries):
		def sample(faceZeros, cellZeros):
			try:
				return np.array(SubLanczosKernelSample(
					self.matrices.coboundary, cellZeros, faceZeros, self.field, maxTries=maxTries
				), dtype=FINT)
			except Exception as e:
				raise NumericalInstabilityWarning(e)
			
		self.sample = sample


	def _initial(self):
		"""
		Computes an initial state for the model's Complex.

		Returns:
			A numpy array of spin assignments (mod p).
		"""
		return self.RNG.integers(
			0, high=self.field, dtype=FINT, size=self.faces
		)
	

	def _proposal(self, time):
		"""
		Proposal scheme for the Nienhuis model.

		Args:
			time (int): Step in the chain.

		Returns:
			A numpy array of spin assignments (mod p).
		"""
		# Choose which cells and faces to include.
		cells = self.RNG.uniform(size=self.cells)
		faces = self.RNG.uniform(size=self.faces)
		threshCells = (cells < self.prob.cell).nonzero()[0]
		threshFaces = (faces < self.prob.face).nonzero()[0]
		q = self.field

		# Evaluate the current spin configuration (cochain) on the boundary
		# of our chosen cells, and find which to include.
		cellBoundary = self.complex.Boundary[self.dimension]
		cellCoefficients = self.spins[cellBoundary]
		cellCoefficients[:,1::2] = -cellCoefficients[:,1::2]%q
		cellSums = cellCoefficients.sum(axis=1)%q
		cellZeros = np.nonzero(cellSums[threshCells] == 0)[0]

		# Decide which faces (by default, edges) to include, take a submatrix
		# of the coboundary matrix by including only the rows corresponding to
		# included cells and columns corresponding to included faces, then
		# sample new spins.
		faceCoefficients = np.intersect1d((self.spins==0).nonzero()[0], threshFaces)
		faceZeros = np.setdiff1d(self.faceIndices, threshFaces)
		subspins = self.sample(faceZeros, cellZeros)

		# Organize spin data.
		spins = self.spins.copy()
		spins[faceZeros] = subspins
		spins[faceCoefficients] = 0

		# Compute energies.
		faceEnergy = (self.faces - np.count_nonzero(spins==0))/self.faces
		cellEnergy = (self.cells - np.count_nonzero(cellSums==0))/self.cells

		return spins, faceEnergy, cellEnergy


	def _assign(self, cocycle):
		"""
		Updates mappings from faces to spins and cubes to occupations.

		Args:
			cocycle (np.array): Cocycle on the subcomplex.
		
		Returns:
			None.
		"""
		self.spins = cocycle
