
# distutils: language=c++

from .common cimport FFINT, FLAT, TABLE, FLATCONTIG, TABLECONTIG

from libcpp.vector cimport vector as Vector

cdef extern from "FastSample.h":
	Vector[int] FastSample(Vector[int], int M, int N, int p);


cpdef Vector[int] Fast(TABLE A, int p) noexcept:
	# Instantiate.
	cdef int i, j, k, t, M, N, L;
	cdef double density;
	cdef Vector[int] _coboundary, coboundary;
	cdef FFINT q;

	M = A.shape[0];
	N = A.shape[1];
	density = 0.1 if M*N > 10e6 else M*N/2
	L = (int)(density*M*N);

	# Put the matrix into a "better" format so we only have to scan over it once.
	t = 0;
	_coboundary = Vector[int](L, 0);

	# Fill the sparse format.
	for i in range(M):
		for j in range(N):
			q = A[i,j];
			if 0 < q and q < p:
				_coboundary[t] = i;
				_coboundary[t+1] = j;
				_coboundary[t+2] = q;
				t += 3;
	
	coboundary = Vector[int](t);
	for k in range(t): coboundary[k] = _coboundary[k];
	
	return FastSample(coboundary, M, N, p);
