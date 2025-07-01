
# distutils: language=c++

from ..common cimport MINT, INDEXFLAT
from .LinBoxMethods cimport LanczosKernelSample as _LanczosKernelSample

from libcpp.vector cimport vector as Vector
from libcpp.unordered_set cimport unordered_set as Set
from libcpp.unordered_map cimport unordered_map as Map

# cdef extern from "FastSample.h":
# 	Vector[int] FastSample(Vector[int], int M, int N, int p);


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

	for i in range(0, coboundary.shape[0], 3):
		cube = coboundary[i];

		if zeroset.contains(cube):
			submatrix[t] = zeromap[cube];
			submatrix[t+1] = coboundary[i+1];
			submatrix[t+2] = coboundary[i+2];
			t += 3

	return _LanczosKernelSample(submatrix, L, columns, field, maxTries);
