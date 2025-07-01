
import numpy as np


def autocorrelation(data):
	"""
	Computes the autocorrelation of a given observable.

	Args:
		data (Iterable): An iterable, indexed by sample times, containing data
			from a given observable.
	
	Returns:
		An `np.ndarray` of autocorrelation data.
	"""

	# Expected value (i.e. sample mean, which converges to the expectation by
	# LLN).
	mu = data.mean();
	normalized = data-mu;
	N = len(data);

	autocorrs = np.array([
		np.dot(normalized[t:], normalized[:N-t])*(1/N) for t in range(N)
	])

	return autocorrs/autocorrs[0]
