
import numpy as np
cimport numpy as np


def totalEnergy(state, model):
	r"""
	Computes the _total energy_ (or, more precisely, the _net energy_) \( \mathcal E \)
	of a given state produced by the given model. The total energy is given by \[ \mathcal E (\sigma)
	= \sum_{x \in X} (-1)^{1-\mathbf 1\left[ \sigma(\partial x) = 0 \right]} =
	\sum_{x \in X} \mathbf 1\left[ \sigma(\partial x) = 0 \right] -
	\sum_{x \in X} \mathbf 1\left[ \sigma(\partial x) \neq 0 \right], \] where
	the \(x \in X\) are \(d\)-cells and \(\sigma\) is a spin configuration on
	\((d-1)\)-cells.

	Args:
		state (tuple): The state from the Markov chain on the given `model`.
		model (Model): The Model (e.g. `ateams.models.SwendsenWang`) from which
			data is collected.

	Returns:
		The total (net) energy.
	"""
	cdef int[:] spins, occupied, satisfied;
	cdef int q;

	spins, occupied, satisfied = state
	q = model.field

	boundary = model.complex.Boundary[model.dimension]
	coefficients = spins[boundary]
	coefficients[:,1::2] = -coefficients[:,1::2]%q
	sums = coefficients.sum(axis=1)%q

	agree = np.nonzero(sums == 0)[0]
	disagree = len(boundary)-agree

	return agree-disagree
	

def occupancy(state, model=None):
	r"""
	Computes the _bond occupancy_ \[ \mathcal N(\sigma, \omega) = \sum_{x \in X}
	\mathbf 1[ \omega(x) = 1 ]\] of the state produced by the given model, where
	\(\omega\) is a bond(/plaquette) configuration. (Just counts the number of
	\(d\)-cells included by independent percolation conditioned on a spin configuration
	\(\sigma\).) Note: at every step of a simulation, each Model returns a
	binary vector over the \(d\)-cells indicating which are in/excluded; the
	occupancy can be computed (as is done here) by `.sum()`ming on this vector.

	Args:
		state (tuple): The state from the Markov chain on the given `model`.
		model (Model=None): The Model (e.g. `ateams.models.SwendsenWang`) from
			which data is collected. Doesn't do anything here.

	Returns:
		The bond occupancy.
	"""
	spins, occupied, satisfied = state
	return occupied.sum()


def connectivity(state, model):
	pass
