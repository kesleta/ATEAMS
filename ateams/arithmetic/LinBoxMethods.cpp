
#include <linbox/linbox-config.h>
#include <linbox/solutions/solve.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/ring/modular.h>

#include <iostream>

#include "LinBoxMethods.h"

using namespace LinBox;
using namespace std;

typedef Givaro::Modular<int> Field;
typedef SparseMatrix<Field, SparseMatrixFormat::SparseSeq> FieldMatrix;
typedef DenseVector<Field> FieldVector;


bool containsNonzero(FieldVector X) {
	// Check whether there are nonzero elements in the VVector.
	bool contains = false;
	int t = 0;

	for (FieldVector::const_iterator it = X.begin(); it != X.end(); ++it) {
		if (*it > 0) {
			contains = true;
			break;
		}
	}

	return contains;
}


FieldMatrix FieldFill(vector<int> coboundary, int M, int N, Field F) {
	// Construct the sparse coboundary matrix.
	FieldMatrix A(F, M, N);

	for (int t = 0; t < coboundary.size(); t += 3) {
		Field::Element q;
		int i, j;

		i = coboundary[t];
		j = coboundary[t+1];
		F.init(q, coboundary[t+2]);
		A.setEntry(i,j,q);
	}

	return A;
}


vector<int> populate(FieldMatrix A, FieldVector X) {
	// Populate a vector with entries from a LinBox vector.
	std::vector<int> x(A.coldim());
	
	for (size_t k = 0; k < A.coldim(); k++) {
		x[k] = X.getEntry(k);
	}

	return x;
}


vector<int> LanczosKernelSample(vector<int> coboundary, int M, int N, int p, int maxTries) {
	// Construct the finite field and construct the matrix.
	Field F(p);
	FieldMatrix A = FieldFill(coboundary, M, N, F);
	FieldVector X(F, A.coldim()), b(F, A.rowdim());

	// Re-sample until we get something other than the zero vector (up to four
	// tries).
	int t = 0;

	while (!containsNonzero(X) and t < maxTries) {
		solve(X, A, b, Method::Lanczos());
		t++;
	}

	return populate(A, X);
}


typedef unordered_set<int> Set;
typedef map<int, Field::Element> Column;
typedef vector<Column> BoundaryMatrix;


BoundaryMatrix FillBoundaryMatrix(vector<int> boundary, Field F, int L) {
	BoundaryMatrix B(L);

	for (int t = 0; t < boundary.size(); t += 3) {
		Field::Element q;
		int i, j;

		i = boundary[t];
		j = boundary[t+1];
		F.init(q, boundary[t+2]);
		B[j][i] = q;
	}

	return B;
}


int dim(int index, vector<int> breaks) {
	int i = 0;

	for (; i < breaks.size()-1; i++) {
		if (index < breaks[i+1]) return i;
	}

	return i;
}

int low(Column column) {
	return column.rbegin()->first;
}


Set ComputePercolationEvents(vector<int> boundary, vector<int> filtration, int p, vector<int> breaks) {
	// Construct the finite field and build the boundary matrix. Similarly to Chen and
	// Kerber (2011), we store each column of the boundary matrix as a balanced
	// binary search tree (as implemented by the C++ standard library). Chen
	// and Kerber (and their successors at PHAT) take advantage of the fact that
	// they are computing coefficients over Z/2Z and need only store the indices
	// of the nonzero entries in each column, because those entries are only ever
	// 1. Here, we have to store the entries as well. For that, we use the standard
	// C++ map, which takes column indices to Givaro::Modular<int> finite field
	// elements for fast arithmetic.
	int cellCount = filtration.size();
	int topDimension = breaks.size()-1;
	int lowIndex, row;

	Field F(p);
	Field::Element q, r, s;
	BoundaryMatrix B = FillBoundaryMatrix(boundary, F, cellCount);
	vector<int> nextColumn = vector<int>(cellCount, 0);
	Column column, replacement, lower;

	// TODO: dimensionality lookup; should probably be done with breaks/tranches
	// but we can figure that out later.
	for (int d = topDimension; d > 0; d--) {
		for (int j = 0; j < cellCount; j++) {
			// If we aren't of the right dimension, don't do anything; I think
			// we can probably improve the efficiency here, since we know when
			// cells of whatever dimension will be placed in the filtration.
			// That should cut down on time quite a bit!
			if (dim(j, breaks) != d) continue;

			// Otherwise, while the current column isn't zero, do some column
			// reduction.
			column = B[j];
			replacement = Column();
			lowIndex = low(column);

			while (!column.empty() & nextColumn[lowIndex] != 0) {
				lower = B[lowIndex];
				q = lower[lower.rend()->first];

				// Row arithmetic. TODO: put this in an actual separate function
				// so it's easier to debug. If this is faster than my old stuff,
				// I'm going to shit myself.
				for (auto it = column.begin(); it != column.end(); ++it) {
					// Hedge that, most of the time, we won't actually be doing
					// that much arithmetic.
					row = it->first;
					r = it->second;

					// If we *do* have a shared element in that row, go nuts.
					if (auto search = lower.find(row); search != lower.end()) {
						s = search->second;
						replacement[row] = r-(1/q)*s;
						lower.erase(row);
					} else {
						replacement[row] = r;
					}
				}
			}
		}
	}

	return Set();
}
