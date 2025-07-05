
import numpy as np

from ..common import MINT, FINT, Matrices
from ..statistics import constant
from .Model import Model


class Glauber(Model):
	name = "Glauber"
	
	def __init__(
			self, C, dimension=1, temperature=constant(-0.6), initial=None
		):
		"""
		Initializes Glauber dynamics on the Potts model.

		Args:
			C: The `Complex` object on which we'll be running experiments.
			dimension (int=1): Dimension on which we'll be building experiments.
			temperature (Callable): A temperature schedule function which
				takes a single positive integer argument `t`, and returns the
				scheduled temperature at time `t`.
			initial (np.ndarray): A vector of spin assignments to components.
		"""
		self.complex = C
		self.dimension = dimension
		self.temperature = temperature

		# Useful values to have later.
		self.cells = len(self.complex.Boundary[self.dimension])
		self.faces = len(self.complex.Boundary[self.dimension-1])
		self.closed = 1

		# Seed the random number generator.
		self.RNG = np.random.default_rng()

		# If no initial spin configuration is passed, initialize.
		if not initial: self.spins = self.initial()
		else: self.spins = (initial%self.complex.field).astype(FINT)
	

	def initial(self):
		"""
		Computes an initial state for the model's Complex.

		Returns:
			A numpy `np.array` representing a vector of spin assignments.
		"""
		return self.RNG.integers(
			0, high=self.complex.field, dtype=MINT, size=self.faces
		)
	

	def proposal(self, time):
		"""
		Proposal scheme for generalized Glauber dynamics on the Potts model:
		uniformly randomly chooses a face in the complex, flips the face's spin,
		and returns the corresponding cocycle.

		Args:
			time (int): Step in the chain.

		Returns:
			A NumPy array representing a vector of spin assignments.
		"""
		# Flip.
		spins = self.spins[::]
		flip = self.RNG.integers(self.faces)
		spins[flip] = (self.RNG.integers(self.complex.field, dtype=FINT))

		boundary = self.complex.Boundary[self.dimension]
		q = self.complex.field

		# Evaluate the current spin assignment (cochain).
		evaluation = spins[boundary]
		evaluation[:, 1::2] = -evaluation[:, 1::2]%q
		sums = evaluation.sum(axis=1)%q

		opens = np.nonzero(sums == 0)[0]
		closed = -(self.faces-opens.shape[0])

		# If we don't know how many "closed" faces there are, compute it.
		if self.closed > 0:
			evaluation = self.spins[boundary]
			evaluation[:,1::2] = -evaluation[:,1::2]%q
			sums = evaluation.sum(axis=1)%q
			
			opens = np.nonzero(sums == 0)[0]
			self.closed = -(self.faces-opens.shape[0])

		# Compute the energy and decide if this is an acceptable transition.
		energy = np.exp(self.temperature(time)*(self.closed-closed))

		if self.RNG.uniform() < min(1, energy):
			self.closed = closed
			self.open = opens
			self.spins = spins

		satisfied = np.zeros(self.cells, dtype=FINT)
		satisfied[self.open] = 1

		return self.spins, satisfied
	

	def assign(self, cochain):
		"""
		Updates mappings from faces to spins and cubes to occupations.
		
		Args:
			cochain (np.array): Cochain on the complex.
		
		Returns:
			None.
		"""
		self.spins = cochain


