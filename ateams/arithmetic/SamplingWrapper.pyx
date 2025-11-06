
# distutils: language=c++

from ..common cimport int, INDEXFLAT, Index, Set, Map, bool
from .Sampling cimport (
	ReducedKernelSample as _ReducedKernelSample
)


cdef Index RowSubmatrix(INDEXFLAT A, INDEXFLAT includeRows, int columns):
	# Keep a set of rows to include.
	cdef int i, t, cube;
	cdef Set zeroset = Set();
	cdef Map zeromap = Map();
	cdef Index submatrix = Index();

	for i in range(includeRows.shape[0]):
		zeroset.insert(includeRows[i]);
		zeromap[includeRows[i]] = i;

	# Create an empty submatrix; iterate over the original matrix, inserting
	# as we go.
	for t in range(0, A.shape[0], 3):
		# Check whether the submatrix includes this cube (which corresponds to
		# the row index).
		cube = A[t];

		if zeroset.contains(cube):
			submatrix.push_back(zeromap[cube]);
			submatrix.push_back(A[t+1]);
			submatrix.push_back(A[t+2]);
	
	return submatrix


cdef Index SubZeroSubmatrix(INDEXFLAT A, INDEXFLAT zeroRows, INDEXFLAT zeroColumns):
	cdef int i, t, L, M, N;
	cdef int cell, face;
	cdef Index submatrix;
	cdef Set zeroRowSet, zeroColumnSet;
	cdef Map zeroRowMap, zeroColumnMap;

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
	submatrix = Index();

	for i in range(0, A.shape[0], 3):
		cell = A[i];
		face = A[i+1]

		if zeroRowSet.contains(cell) and zeroColumnSet.contains(face):
			submatrix.push_back(zeroRowMap[cell]);
			submatrix.push_back(zeroColumnMap[face]);
			submatrix.push_back(A[i+2]);

	return submatrix


cpdef Index SubReducedKernelSample(
		INDEXFLAT coboundary, INDEXFLAT zeroRows, INDEXFLAT zeroColumns,
		int p, bool verbose
	) except *:
	"""
	Uses the SpaSM sparse Gaussian elimination strategy to sample uniformly from
	the kernel of the submatrix of the coboundary matrix.

	Args:
		coboundary: Memory-contiguous one-dimensional coboundary array,
			like that of the `Cubical` Complex object.
		zeros: Memory-contiguous array of row indices that will be included
			in the coboundary submatrix.
		p: Characteristic of the field \(\mathbb Z/p\mathbb Z\).

	Returns:
		A C++ `std::vector` with spin assignments.
	"""
	cdef Index submatrix = SubZeroSubmatrix(coboundary, zeroRows, zeroColumns);
	return _ReducedKernelSample(submatrix, zeroRows.shape[0], zeroColumns.shape[0], p, verbose);


cpdef Index ReducedKernelSample(
		INDEXFLAT coboundary, INDEXFLAT zeros, int faces, int columns, int p,
		bool verbose
	) except *:
	"""
	Uses the SpaSM sparse Gaussian elimination strategy to sample uniformly from
	the kernel of the submatrix of the coboundary matrix.

	Args:
		coboundary: Memory-contiguous one-dimensional coboundary array,
			like that of the `Cubical` Complex object.
		zeros: Memory-contiguous array of row indices that will be included
			in the coboundary submatrix.
		faces: The number of faces per plaquette.
		columns: Integer representing the number of columns in the coboundary
			submatrix; should be the total number of faces in the complex.
		p: Characteristic of the field \(\mathbb Z/p\mathbb Z\).

	Returns:
		A C++ `std::vector` with spin assignments.
	"""
	cdef Index submatrix = RowSubmatrix(coboundary, zeros, columns);
	return _ReducedKernelSample(submatrix, zeros.shape[0], columns, p, verbose);
