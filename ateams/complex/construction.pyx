
# distutils: language=c++

from libcpp.vector cimport vector as Vector

from ..common cimport INDEX, MINT


cpdef Vector[Vector[MINT]] boundaryMatrices(MINT[:,:] cubes, MINT[:] coefficients):
	cdef INDEX boundary;
	cdef int i, j, face, M, N;
	cdef Vector[MINT] Boundary, Coboundary;
	cdef Vector[Vector[MINT]] Matrices;

	# Instantiate vectors so we don't have to allocate memory all at once.
	Boundary = Vector[MINT]();
	Coboundary = Vector[MINT]();
	Matrices = Vector[Vector[MINT]](2);

	M = cubes.shape[0];
	N = cubes.shape[1];

	# Build the boundary and coboundary matrices.
	for i in range(M):
		boundary = cubes[i];

		for j in range(N):
			face = boundary[j];
			Boundary.push_back(face);
			Boundary.push_back(i);
			Boundary.push_back(coefficients[j]);

			Coboundary.push_back(i);
			Coboundary.push_back(face);
			Coboundary.push_back(coefficients[j]);

	# Return a 2D array: the first row is the boundary matrix, the second is
	# the coboundary matrix (both flattened).
	Matrices[0] = Boundary;
	Matrices[1] = Coboundary;

	return Matrices;


cpdef Vector[int] fullBoundaryMatrix(list boundary, dict coefficients):
	cdef int M, i, N, j, face;
	cdef Vector[int] Boundary = Vector[int]();

	M = len(boundary);
	
	for i in range(M):
		faces = boundary[i];

		# If we find a vertex, we just push back the required entries (a 1 at the
		# (i,i)th entry) and be done.
		try: N = len(faces);
		except:
			Boundary.push_back(i);
			Boundary.push_back(i);
			Boundary.push_back(0);
			continue;

		# Otherwise, we iterate over the faces, adding the required entries.
		for j in range(N):
			face = faces[j];
			Boundary.push_back(face);
			Boundary.push_back(i);
			Boundary.push_back(<int>coefficients[N//2][j]);

	return Boundary;
