
#include <linbox/linbox-config.h>
#include <linbox/solutions/solve.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/ring/modular.h>

#include <iostream>
#include <chrono>
#include <thread>

#include "LinBoxMethods.h"

using namespace std;

typedef vector<int> Vector;
typedef Givaro::Modular<int> Field;
typedef LinBox::SparseMatrix<Field, LinBox::SparseMatrixFormat::SparseSeq> FieldMatrix;
typedef LinBox::DenseVector<Field> FieldVector;


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


FieldMatrix FieldFill(Vector coboundary, int M, int N, Field F) {
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


Vector populate(FieldMatrix A, FieldVector X) {
	// Populate a vector with entries from a LinBox vector.
	Vector x(A.coldim());
	
	for (size_t k = 0; k < A.coldim(); k++) {
		x[k] = X.getEntry(k);
	}

	return x;
}


Vector LanczosKernelSample(Vector coboundary, int M, int N, int p, int maxTries) {
	// Construct the finite field and construct the matrix.
	Field F(p);
	FieldMatrix A = FieldFill(coboundary, M, N, F);
	FieldVector X(F, A.coldim()), b(F, A.rowdim());

	// Re-sample until we get something other than the zero vector (up to four
	// tries).
	int t = 0;
	size_t rank;

	while (!containsNonzero(X) and t < maxTries) {
		LinBox::Method::Lanczos LANK;
		LANK.preconditioner = LinBox::Preconditioner::FullDiagonal;
		solve(X, A, b, LANK);
		t++;
	}

	return populate(A, X);
}


typedef set<int> Set;
typedef map<int, Field::Element> Column;
typedef map<int,int> Map;
typedef vector<Column> BoundaryMatrix;


BoundaryMatrix FillBoundaryMatrix(Vector boundary, Field F, int L) {
	// Fills a sparse representation of a boundary matrix. Uses standard maps,
	// since they're red-black trees that allow for fast min/max lookups (and
	// are self-balancing on removals); take integer indices to Givaro modular
	// integers.
	BoundaryMatrix B(L);

	for (int t = 0; t < boundary.size(); t += 3) {
		Field::Element q;
		int i, j;

		i = boundary[t];
		j = boundary[t+1];
		F.init(q, boundary[t+2]);
		
		if (q > 0) B[j][i] = q;
	}

	return B;
}


int dim(int index, Vector breaks) {
	// Looks up the dimension of a cell based on the breaks.
	int i = 0;

	for (; i < breaks.size()-1; i++) {
		if (index < breaks[i+1]) return i;
	}

	return i;
}


int youngestOf(Column column) {
	// Gets the "youngest" (largest-indexed) cell in the column.
	return column.rbegin()->first;
}


void printMap(Map column) {
	cout << "{" << endl;
	for (auto it = column.begin(); it != column.end(); ++it) {
		cout << "\t" << it->first << ": " << it->second << endl;
	}
	cout << "}" << endl;
}


void printColumn(Column column) {
	cout << "{" << endl;
	for (auto it = column.begin(); it != column.end(); ++it) {
		cout << "\t" << it->first << ": " << it->second << endl;
	}
	cout << "}" << endl;
}

void printSet(Set S) {
	cout << "{";
	for (auto it = S.begin(); it != S.end(); ++it) {
		cout << *it << ", ";
	}
	cout << "}" << endl;
}


void printBoundary(Vector boundary) {
	cout << "[" << endl;
	for (int i=0; i < boundary.size(); i+=3) {
		cout << "  [ " << boundary[i] << " " << boundary[i+1] << " " << boundary[i+2] << " ]" << endl;
	}
	cout << "]" << endl;
}


Set ComputePercolationEvents(
		Vector boundary, Vector filtration, int homology, int p,
		Vector breaks
	) {
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
	int topDimension = (homology < breaks.size() ? homology+1 : breaks.size());
	int row, j, high;

	Field F(p);
	BoundaryMatrix Boundary = FillBoundaryMatrix(boundary, F, cellCount);
	Vector nextColumn = Vector(cellCount, 0);
	Column column, younger;
	Set erase, skip, marked;

	// Initialize some things.
	// TODO: PRE-MARK THINGS WE KNOW WILL BE MARKED EVENTUALLY
	Field::Element q, r, s, inv, prod, result;
	F.init(inv);
	F.init(prod);
	F.init(result);

	// TODO: dimensionality lookup; should probably be done with breaks/tranches
	// but we can figure that out later.
	for (int d = topDimension; d > homology-1; d--) {
		// Since we know the cells will always be added in order of dimension,
		// we know specifically which range to search (saving us a lot of time).
		high = (d+2 >= breaks.size() ? cellCount : breaks[d+2]);
		for (int j = breaks[d-1]; j < high; j++) {
			// If we aren't of the right dimension, don't do anything; I think
			// we can probably improve the efficiency here, since we know when
			// cells of whatever dimension will be placed in the filtration.
			// That should cut down on time quite a bit!
			if (dim(j, breaks) != d) continue;

			// Otherwise, while the current column isn't zero, do some column
			// reduction.
			column = Boundary[j];

			while (!column.empty() && nextColumn[youngestOf(column)] != 0) {
				younger = Boundary[nextColumn[youngestOf(column)]];
				q = younger[youngestOf(column)];

				erase = Set();
				skip = Set();

				// Row arithmetic. TODO: put this in an actual separate function
				// so it's easier to debug. If this is faster than my old stuff,
				// I'm going to shit myself.
				for (auto it = column.begin(); it != column.end(); ++it) {
					// Hedge that, most of the time, we won't actually be doing
					// that much arithmetic.
					row = it->first;
					r = it->second;

					// If we share a row, compute the result; if the result is
					// nonzero, overwrite whatever's in there already.
					if (younger.count(row) > 0) {
						s = younger[row];
						
						F.inv(inv, q);
						F.mul(prod, inv, s);
						F.sub(result, r, prod);
						column[row] = result;

						if (result == 0) erase.insert(row);

						skip.insert(row);
					}
				}

				// Go through and erase the rows that were zeroed out.
				for (auto it = erase.begin(); it != erase.end(); ++it) column.erase(*it);

				// Go through and eliminate the remaining rows, skipping ones
				// we've already computed.
				for (auto it = younger.begin(); it != younger.end(); ++it) {
					row = it->first;
					s = it->second;

					if (skip.count(row) > 0) continue;

					F.inv(inv, q);
					F.mul(prod, inv, s);
					F.neg(result, prod);
					column[row] = result;
				}
			}

			// If there is no younger column to add and the column is nonempty,
			// then we've found a pivot column. On the other hand, if there's no
			// younger column to add and the column is empty, then this column
			// (corresponding to the boundary of the jth cell) is in the basis
			// of the kernel --- congrats, it's a cycle!
			if (!column.empty()) {
				nextColumn[youngestOf(column)] = j;
				Boundary[youngestOf(column)] = Column();
			} else {
				marked.insert(j);
			}
		}
	}

	// Find the essential birth times by checking whether the column is marked
	// (i.e. is a cycle) and has no younger columns to add. (For some reason,
	// 0 gets left out here. Not sure why...)
	Set essential = Set();
	essential.insert(0);
	
	for (auto it = marked.begin(); it != marked.end(); it++) {
		if (nextColumn[*it] == 0) essential.insert(*it);
	}

	return essential;
}
