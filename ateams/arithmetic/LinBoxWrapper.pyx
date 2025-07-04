
# distutils: language=c++

from ..common cimport MINT, INDEXFLAT
from .LinBoxMethods cimport (
	LanczosKernelSample as _LanczosKernelSample,
	ComputePercolationEvents as _ComputePercolationEvents
)

from libcpp.vector cimport vector as Vector
from libcpp.set cimport set as Set
from libcpp.unordered_map cimport unordered_map as Map


cdef Vector[int] Vectorize(INDEXFLAT A) noexcept:
	cdef int i, L;
	cdef Vector[int] B;

	L = A.shape[0];
	B = Vector[int](L);

	for i in range(L): B[i] = A[i];

	return B;


cdef Vector[int] ZeroSubmatrix(INDEXFLAT A, INDEXFLAT zeros, int faces, int columns):
	cdef int i, t, L, M, N;
	cdef int cube;
	cdef Vector[int] submatrix;
	cdef Set[MINT] zeroset;
	cdef Map[MINT, MINT] zeromap;

	# Initialize zeroset and zeromap.
	L = zeros.shape[0];

	for i in range(L):
		zeroset.insert(zeros[i]);
		zeromap[zeros[i]] = i;

	# Create the submatrix: we know how many cells (rows) are going to be included.
	# Each row has `faces` entries, and each entry requires three cells of the
	# array: two for the indices and one for the value. Since the number of entries
	# is faces*cells, the total number of entries is faces*cells*3.
	M = 3*faces*L;
	t = 0;
	submatrix = Vector[int](M);

	for i in range(0, A.shape[0], 3):
		cube = A[i];

		if zeroset.contains(cube):
			submatrix[t] = zeromap[cube];
			submatrix[t+1] = A[i+1];
			submatrix[t+2] = A[i+2];
			t += 3

	return submatrix


cpdef Vector[int] LanczosKernelSample(
		INDEXFLAT coboundary, INDEXFLAT zeros, int faces, int columns, int field,
		int maxTries=16
	) noexcept:
	"""
	Uses the LinBox (Block?) Lanczos black-box solver to get a uniform
	random element of the kernel of the coboundary.

	Args:
		coboundary: Memory-contiguous one-dimensional coboundary array,
			like that of the `Cubical` Complex object.
		zeros: Memory-contiguous array of row indices that will be included
			in the coboundary submatrix.
		faces: The number of faces per plaquette.
		columns: Integer representing the number of columns in the coboundary
			submatrix; should be the total number of faces in the complex.
		field: Integer representing the characteristic of the finite field
			over which computations are done.
		maxTries: How many times do we try to get a nonzero solution from
			the black-box solver? Default is 16.

	Returns:
		A C++ `std::vector` with spin assignments.
	"""
	cdef Vector[int] submatrix = ZeroSubmatrix(coboundary, zeros, faces, columns);
	return _LanczosKernelSample(submatrix, zeros.shape[0], columns, field, maxTries);


cpdef Set[int] ComputePercolationEvents(INDEXFLAT boundary, INDEXFLAT filtration, int homology, int p, INDEXFLAT breaks) noexcept:
	"""
	Uses a variant of the Chen/Kerber (2011) and PHAT (2017) twist_reduce algorithm
	to compute the persistent homology of the complex specified by the flat boundary
	matrix and the filtration, over the field Z/pZ.

	Args:
		boundary: Full boundary matrix (i.e. `Cubical.matrices.full`).
		filtration: A list of indices specifying the order in which the cells
			in the complex will be added. The assumed order is adding cells
			in order of dimension.
		homology: The homology group for which persistence is computed.
		p: Characteristic of the field \(\mathbb{Z}/p\mathbb{Z}\).
		breaks: Index locations for cells of each degree (i.e. `Cubical.breaks`).
	"""
	cdef Vector[int] _boundary = Vectorize(boundary);
	cdef Vector[int] _filtration = Vectorize(filtration);
	cdef Vector[int] _breaks = Vectorize(breaks);

	return _ComputePercolationEvents(_boundary, _filtration, homology, p, _breaks);
