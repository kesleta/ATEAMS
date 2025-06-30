 
import numpy as np

from ..arithmetic import evaluateCochain, FINT, KernelSample, MatrixReduction, Fast
from ..stats import constant
from .Model import Model


class SwendsenWang(Model):
	name = "SwendsenWang"
	
	def __init__(
			self, L, temperatureFunction=constant(-0.6), initial=None, sparse=True,
			parallel=False, minBlockSize=32, maxBlockSize=64, cores=4, LinBox=False
		):
		"""
		Initializes Swendsen-Wang evolution on the Potts model.

		Args:
			L: The `Lattice` object on which we'll be running experiments.
			temperatureFunction (Callable): A temperature schedule function which
				takes a single positive integer argument `t`, and returns the
				scheduled temperature at time `t`.
			initial (galois.FieldArray): A vector of spin assignments to components.
			sparse (bool=True): Should matrices be formatted sparsely? (Uses C/C++).
			parallel (bool=False): Should matrix computations be done in parallel? (Uses C/C++).
			minBlockSize (int=32): If `parallel` is truthy, this is the smallest
				number of columns processed in parallel.
			maxBlockSize (int=64): If `parallel` is truthy, this is the largest
				number of columns processed in parallel.
			cores (int=4): Number of available CPUs/cores/threads on the machine.
			LinBox (bool=False): Uses fast LinBox routines.
		"""
		self.lattice = L
		self.temperatureFunction = temperatureFunction
		self.faceCells = len(self.lattice.boundary[self.lattice.dimension-1])
		self.cubeCells = len(self.lattice.boundary[self.lattice.dimension])

		self.sparse = sparse
		self.parallel = parallel
		self.minBlockSize = minBlockSize
		self.maxBlockSize = maxBlockSize
		self.cores = cores

		self.coboundary = self.lattice.matrices.coboundary.astype(FINT)
		self.Reducer = MatrixReduction(self.lattice.field.characteristic, parallel, minBlockSize, maxBlockSize, cores)
		self.identity = np.identity(self.coboundary.shape[1], dtype=FINT)
		self.SampleFromKernel = KernelSample
		self.LinBox = LinBox

		# SW defaults.
		self.spins = initial if initial else self.initial()


	def initial(self):
		"""
		Computes an initial state for the model's Lattice.

		Returns:
			A Galois `Array` representing a vector of spin assignments.
		"""
		return self.lattice.field.Random(self.faceCells)
	

	def proposal(self, time):
		"""
		Proposal scheme for generalized Swendsen-Wang evolution on the Potts model.

		Args:
			time (int): Step in the chain.

		Returns:
			A Galois `Array` representing a vector of spin assignments.
		"""
		# Compute the probability of choosing any individual cube in the complex.
		self.temperature = self.temperatureFunction(time)
		p = 1-np.exp(self.temperature)
		assert 0 <= p <= 1

		# Choose cubes to include; in effect, this just does a boatload of indexing.
		uniforms = np.random.uniform(size=self.cubeCells)
		include = (uniforms < p).nonzero()[0]
		boundary = self.lattice.boundary[self.lattice.dimension][include]
		boundaryValues = evaluateCochain(boundary, self.spins)
		zeros = (boundaryValues == 0).nonzero()[0]

		satisfied = np.zeros(len(self.lattice.boundary[self.lattice.dimension])).astype(int)
		satisfied[zeros] = 1

		# Uniformly randomly sample a cocycle on the sublattice admitted by the
		# chosen edges; reconstruct the labeling on the entire lattice by
		# subbing in the values of c which differ from existing ones.
		if not self.LinBox:
			spins = self.SampleFromKernel(self.Reducer, self.coboundary.take(zeros, axis=0))
		else:
			spins = Fast(self.coboundary.take(zeros, axis=0), self.lattice.field.characteristic)

		spins = self.lattice.field(spins)

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
