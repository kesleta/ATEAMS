
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
	mu = data.mean()
	normalized = data-mu
	M = data.shape[0]

	autocorrs = np.empty(M)
	for t in range(M): autocorrs[t] = np.dot(normalized[t:], normalized[:M-t])*(1/M)

	return autocorrs/autocorrs[0]

