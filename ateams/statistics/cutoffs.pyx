

cpdef int rectangular(Z, c=7):
	r"""
	Computes the rectangular cutoff for the integrated autocorrelation time estimator
	\(Z\) --- this cutoff is (approximately) the number of "good" estimates of the
	integrated autocorrelation time before random walk-y behavior occurs, as in
	Ossola and Sokal (2004).

	Args:
		Z (np.array): NumPy array of integrated autocorrelation values.
		c (int=7): Appropriate scaling constant.
	"""
	cdef int M = 0;

	while M < c*Z[M]: M += 1

	return M
