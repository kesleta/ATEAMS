
import numpy as np
cimport numpy as np


def integrated(X, int M=-1):
	r"""
	Computes (an estimate of) the _integrated autocorrelation time_ \(\tau_X\) of the observable \(X\)
	by \[ \tau_X = \frac 12 \sum_{t=-M}^M \overline \rho_X(t), \] where \(\overline \rho_X(t)\)
	is the normalized autocorrelation for \(X\).

	Args:
		X (np.array): A NumPy array of numerical values.
		M (double=-1): The autocorrelation window width. Defaults to \(M = |X|\),
			which uses all observed data.

	Returns:
		(An estimate of) the integrated autocorrelation time.
	"""
	if M < 0: M = X.shape[0]

	auto = normalized(X)
	return auto[:M].sum()


def unnormalized(X):
	r"""
	Computes the vector \(\rho_X\) of unnormalized autocorrelations \(\rho_X(t) \)
	for the set of (numerical) observations \(X\) by
	\[ \rho_X(t) = \frac{1}{N-|t|} \sum_{i=1}^{N-|t|}(X_i - \overline X)(X_{i+t} - \overline X), \]
	where \(\overline X\) is the sample mean and \(N = |X|\) is the number of
	samples in the observed data.

	Args:
		X (np.array): A NumPy array of numerical values.

	Returns:
		A NumPy array of unnormalized autocorrelation values.
	"""
	cdef double mu = X.mean();
	cdef int t, N = X.shape[0]
	adj = X-mu
	var = np.empty(N)

	for t in range(N):
		u = abs(t)
		var[t] = (1/(N-u))*(adj[:N-u]*adj[u:N]).sum()

	return var


def normalized(X):
	r"""
	Computes the vector \(\overline \rho_X\) of normalized autocorrelations
	\(\overline \rho_X(t) = \rho_X(t)/\rho_X(0) \) for the observable \(X\).

	Args:
		X (np.array): A NumPy array of numerical values.

	Returns:
		A NumPy array of unnormalized autocorrelation values.
	"""
	unn = unnormalized(X)
	return unn/unn[0]
