
# distutils: language=c++

from ..common cimport int, INDEXFLAT, Vector, Set, AnyMap
from .LinBoxMethods cimport (
	LanczosKernelSample as _LanczosKernelSample
)


cdef Vector ZeroSubmatrix(INDEXFLAT A, INDEXFLAT zeros, int faces, int columns):
	cdef int i, t, L, M, N;
	cdef int cube;
	cdef Vector submatrix;
	cdef Set zeroset;
	cdef AnyMap[int, int] zeromap;

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
	submatrix = Vector(M);

	for i in range(0, A.shape[0], 3):
		cube = A[i];

		if zeroset.contains(cube):
			submatrix[t] = zeromap[cube];
			submatrix[t+1] = A[i+1];
			submatrix[t+2] = A[i+2];
			t += 3

	return submatrix


cdef Vector SubZeroSubmatrix(INDEXFLAT A, INDEXFLAT zeroRows, INDEXFLAT zeroColumns):
	cdef int i, t, L, M, N;
	cdef int cell, face;
	cdef Vector submatrix;
	cdef Set zeroRowSet, zeroColumnSet;
	cdef AnyMap[int, int] zeroRowMap, zeroColumnMap;

	# Initialize the set/map for zero rows.
	L = zeroRows.shape[0];

	for i in range(L):
		zeroRowSet.insert(zeroRows[i]);
		zeroRowMap[zeroRows[i]] = i;
	
	# Do the same for columns.
	L = zeroColumns.shape[0];

	for i in range(L):
		zeroColumnSet.insert(zeroColumns[i]);
		zeroColumnMap[zeroColumns[i]] = i;

	# Create the submatrix. Maybe could do some better memory allocation on this
	# but that's a chore for later Anthony.
	submatrix = Vector();

	for i in range(0, A.shape[0], 3):
		cell = A[i];
		face = A[i+1]

		if zeroRowSet.contains(cell) and zeroColumnSet.contains(face):
			submatrix.push_back(zeroRowMap[cell]);
			submatrix.push_back(zeroColumnMap[face]);
			submatrix.push_back(A[i+2]);

	return submatrix


cpdef Vector SubLanczosKernelSample(
		INDEXFLAT coboundary, INDEXFLAT zeroRows, INDEXFLAT zeroColumns,
		int p, int maxTries=8
	) except *:
	"""
	Uses the LinBox (Block?) Lanczos black-box solver to get a uniform
	random element of the kernel of the coboundary.

	Args:
		coboundary: Memory-contiguous one-dimensional coboundary array,
			like that of the `Cubical` Complex object.
		zeros: Memory-contiguous array of row indices that will be included
			in the coboundary submatrix.
		p: Characteristic of the field \(\mathbb Z/p\mathbb Z\).
		maxTries: How many times do we try to get a nonzero solution from
			the black-box solver? Default is 8.

	Returns:
		A C++ `std::vector` with spin assignments.
	"""
	cdef Vector submatrix = SubZeroSubmatrix(coboundary, zeroRows, zeroColumns);
	return _LanczosKernelSample(submatrix, zeroRows.shape[0], zeroColumns.shape[0], p, maxTries);


cpdef Vector LanczosKernelSample(
		INDEXFLAT coboundary, INDEXFLAT zeros, int faces, int columns, int p,
		int maxTries=8
	) except *:
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
		p: Characteristic of the field \(\mathbb Z/p\mathbb Z\).
		maxTries: How many times do we try to get a nonzero solution from
			the black-box solver? Default is 8.

	Returns:
		A C++ `std::vector` with spin assignments.
	"""
	cdef Vector sample, submatrix = ZeroSubmatrix(coboundary, zeros, faces, columns);
	return _LanczosKernelSample(submatrix, zeros.shape[0], columns, p, maxTries);


# cpdef Set ComputePercolationEvents(INDEXFLAT boundary, INDEXFLAT filtration, int homology, int p, INDEXFLAT breaks) noexcept:
# 	"""
# 	Uses a variant of the Chen/Kerber (2011) and PHAT (2017) twist_reduce algorithm
# 	to compute the persistent homology of the complex specified by the flat boundary
# 	matrix and the filtration over the field \(\mathbb{Z}/p\mathbb{Z}\).

# 	Args:
# 		boundary: Full boundary matrix (i.e. `Cubical.matrices.full`).
# 		filtration: A list of indices specifying the order in which the cells
# 			in the complex will be added. The assumed order is adding cells
# 			in order of dimension.
# 		homology: The homology group for which persistence is computed.
# 		p: Characteristic of the field \(\mathbb{Z}/p\mathbb{Z}\).
# 		breaks: Index locations for cells of each degree (i.e. `Cubical.breaks`).
# 	"""
# 	cdef Vector _boundary = Vectorize(boundary);
# 	cdef Vector _filtration = Vectorize(filtration);
# 	cdef Vector _breaks = Vectorize(breaks);

# 	return _ComputePercolationEvents(_boundary, _filtration, homology, p, _breaks);
