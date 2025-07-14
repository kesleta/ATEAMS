
# distutils: language=c++

from ..common cimport INDEXFLAT, Vectorize, Pairs, Vector
from .PHATMethods cimport PHATComputePersistencePairs as _PHATComputePersistencePairs


cpdef Pairs PHATComputePersistencePairs(INDEXFLAT boundary, INDEXFLAT filtration, int homology, INDEXFLAT breaks) noexcept:
	"""
	Computes the persistence pairs of the complex corresponding to the provided
	boundary matrix and filtration.

	Args:
		boundary (np.array): Flattened boundary matrix given by `Complex.matrices.full`.
		filtration (np.array): Permutation on the order of the columns of the
			boundary matrix. **Here, we assume that only cells of dimension `homology`
			are being permuted. Shuffling the order of cells of different
			dimensions will result in incorrect computations.**
		homology (int): Homology group we're interested in; corresponds to the
			dimension of permuted cells.
		breaks (np.array): Index ranges for cells by dimension, given by
			`Complex.breaks`.

	Returns:
		A list of [birth, death] pairs.
	"""
	cdef Vector _boundary, _filtration, _breaks;
	_boundary = Vectorize(boundary);
	_filtration = Vectorize(filtration);
	_breaks = Vectorize(breaks);

	return _PHATComputePersistencePairs(_boundary, _filtration, homology, _breaks);

