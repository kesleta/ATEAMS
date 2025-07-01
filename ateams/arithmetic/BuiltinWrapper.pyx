
# distutils: language=c++

import numpy as np
cimport numpy as np
from ..common cimport FFINT, FLAT, TABLE, FLATCONTIG, TABLECONTIG
from ..common import FINT
from .MatrixReduction cimport MatrixReduction

from libcpp cimport bool
from libc.stdlib cimport rand, srand
from libc.time cimport time

# Seed the RNG.
srand(time(NULL));


cdef int Random(int MAX) except -1:
	return (rand() % MAX);


cdef FLATCONTIG LinearCombinationReduced(
		MatrixReduction MATRIX,
		TABLECONTIG basis
	) noexcept:
	cdef int i, j, N, M;
	cdef FLATCONTIG result, store;

	M = basis.shape[0];
	N = basis.shape[1];

	result = np.zeros_like(basis[0]);
	store = np.zeros(M, dtype=FINT);

	# Pick (uniform random) coefficients in the field.
	for i in range(M):
		store[i] = MATRIX.negation[Random(MATRIX.characteristic)]

	for j in range(N):
		for i in range(M):
			mult = MATRIX.multiplication[store[i], basis[i,j]]
			result[j] = MATRIX.addition[result[j], mult]

	return result


cpdef TABLECONTIG Kernel(
		MatrixReduction M,
		TABLE A
	) noexcept:
	"""
	Returns a basis for the kernel of `A`.

	Args:
		M (MatrixReduction): A `MatrixReduction` object.
		A (np.array): A 2D `NumPy` array to be reduced.

	Returns:
		A 2D `NumPy` array representing the reduced basis of \(\mathrm{ker}(A)\).
	"""
	cdef TABLECONTIG B, I, inversion, reduced, superreduced, eliminated;
	cdef int m, n, minzero;

	m = A.shape[0];
	n = A.shape[1];

	I = np.identity(n, dtype=FINT)
	B = np.concatenate([A.T, I], axis=1, dtype=FINT);
	M.RREF(B, m);

	minzero = M.HighestZeroRow(m);

	if minzero < 0:
		eliminated = np.empty((2,0), dtype=FINT)
	else:
		inversion = M.ToArray();
		reduced = inversion[:,m:];
		superreduced = reduced[minzero:]
		M.RREF(superreduced)
		eliminated = M.ToArray()

	return eliminated


cdef TABLECONTIG __Kernel(
		TABLECONTIG empty,
		MatrixReduction M,
		TABLE A
	) noexcept:
	cdef TABLECONTIG B, I, inversion, reduced, superreduced, eliminated;
	cdef int m, n, minzero;

	m = A.shape[0];
	n = A.shape[1];

	I = np.identity(n, dtype=FINT)
	B = np.concatenate([A.T, I], axis=1, dtype=FINT);
	M.RREF(B, m);

	minzero = M.HighestZeroRow(m);

	if minzero < 0:
		eliminated = empty;
	else:
		inversion = M.ToArray();
		reduced = inversion[:,m:];
		superreduced = reduced[minzero:]
		eliminated = superreduced

	return eliminated


cpdef np.ndarray[FFINT, ndim=1, negative_indices=False, mode="c"] KernelSample(
		MatrixReduction M, TABLE A
	) noexcept:
	"""
	Draws a uniform random sample from the kernel of `A` using the `MatrixReduction`
	object `M`.

	Args:
		M (MatrixReduction): Contains the RREF form of a matrix.
		A (np.ndarray): 2D `NumPy` array whose kernel is computed and sampled
			from.
	Returns:
		A vector sampled at uniform random from the kernel of the matrix stored in
		`M`.
	"""
	cdef TABLECONTIG empty = np.empty((2,0), dtype=FINT);
	cdef TABLECONTIG basis = __Kernel(empty, M, A);
	return np.asarray(LinearCombinationReduced(M, basis));

