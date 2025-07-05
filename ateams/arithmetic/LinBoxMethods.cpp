
#include <linbox/linbox-config.h>
#include <linbox/solutions/solve.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/ring/modular.h>

#include <iostream>

#include "LinBoxMethods.h"

using namespace std;

typedef vector<int> Vector;
typedef Givaro::Modular<int> Field;
typedef LinBox::SparseMatrix<Field, LinBox::SparseMatrixFormat::SparseSeq> FieldMatrix;
typedef LinBox::DenseVector<Field> FieldVector;


bool containsNonzero(FieldVector X) {
	// Check whether there are nonzero elements in the VVector.

	for (auto it = X.begin(); it != X.end(); ++it) {
		if (*it > 0) { return true; }
	}

	return false;
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

	// Use the sparse Lanczos solver; if there are more rows than columns, we can
	// precondition the matrix.
	LinBox::Method::Lanczos LANC;
	LANC.preconditioner = LinBox::Preconditioner::FullDiagonal;

	// Re-sample until we get something other than the zero vector (up to four
	// tries).
	int t = 0;

	while (!containsNonzero(X) and t < maxTries) {
		solve(X, A, b, LANC);
		t++;
	}

	return populate(A, X);
}


typedef set<int> Set;
typedef map<int, Field::Element> Column;
typedef map<int,int> Map;
typedef map<int,Column> ColumnMap;
typedef vector<Column> BoundaryMatrix;
typedef vector<Set> Mod2BoundaryMatrix;

Mod2BoundaryMatrix Mod2FillBoundaryMatrix(Vector boundary, int L) {
	// Fills a sparse representation of a boundary matrix. Uses standard maps,
	// since they're red-black trees that allow for fast min/max lookups (and
	// are self-balancing on removals); take integer indices to Givaro modular
	// integers.
	Mod2BoundaryMatrix B(L);
	int row, column, q;

	for (int t = 0; t < boundary.size(); t += 3) {
		row = boundary[t];
		column = boundary[t+1];
		q = boundary[t+2];

		// Just a check to exclude vertices.
		if (q != 0) { B[column].insert(row); }
	}

	return B;
}


BoundaryMatrix FillBoundaryMatrix(Vector boundary, Field F, int L) {
	// Fills a sparse representation of a boundary matrix. Uses standard maps,
	// since they're red-black trees that allow for fast min/max lookups (and
	// are self-balancing on removals); take integer indices to Givaro modular
	// integers.
	BoundaryMatrix B(L);
	Field::Element q;
	int row, column;

	for (int t = 0; t < boundary.size(); t += 3) {

		row = boundary[t];
		column = boundary[t+1];
		F.init(q, boundary[t+2]);
		
		if (q > 0) { B[column][row] = q; }
	}

	return B;
}


int Mod2youngestOf(Set column) {
	// Gets the "youngest" (largest-indexed) cell in the column.
	return *column.rbegin();
}


template <typename BalancedStorage>
int youngestOf(BalancedStorage column) {
	// Gets the "youngest" (largest-indexed) cell in the column.
	return column.rbegin()->first;
}


Vector ReindexBoundaryMatrix(Vector &boundary, Vector filtration, int homology, Vector breaks) {
	// Decide when we should be editing rows or columns.
	int cellCount = filtration.size();
	int low = breaks[homology];
	int high = breaks[homology+1];
	int higher = (homology+2 >= breaks.size() ? cellCount : breaks[homology+2]);
	int row, column;

	// Reindex the boundary matrix in-place so we don't have to do it in Python.
	// Eugh.
	Map IndexMap = Map();

	for (int t=low; t < higher; t++) {
		IndexMap[filtration[t]] = t;
	}

	for (int i=0; i < boundary.size(); i+=3) {
		row = boundary[i];
		column = boundary[i+1];

		if (low <= column && column < high) {
			// In this situation, we're editing the *column* --- for example,
			// if the (usual) 9th element was placed 10th, then we have to change
			// the 9th column to the 10th one.
			boundary[i+1] = IndexMap[column];
		} else if (high <= column && column < higher) {
			// Here, we're editing the *row* --- since we're in a higher-dimensional
			// cell, we need to point toward the correct (re-indexed) row.
			boundary[i] = IndexMap[row];
		}
	}
	
	return boundary;
}


Set Mod2ComputePercolationEvents(
		Vector boundary, Vector filtration, int homology, Vector breaks
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

	Vector _boundary = ReindexBoundaryMatrix(boundary, filtration, homology, breaks);
	Mod2BoundaryMatrix Boundary = Mod2FillBoundaryMatrix(_boundary, cellCount);

	Vector nextColumnAdded = Vector(cellCount, 0);
	Set cell, youngest, sym, marked;

	int face, high, numBreaks = breaks.size();
	marked = Set();

	for (int d = topDimension; d > homology-1; d--) {

		high = (d+1 >= numBreaks ? cellCount : breaks[d+1]);

		for (int j = breaks[d]; j < high; j++) {
			// If we're of the wrong dimension, keep going.
			// if (dim(j, breaks) != d) { continue; }
			cell = Boundary[j];

			while (!cell.empty() && nextColumnAdded[Mod2youngestOf(cell)] != 0) {
				// Get the "youngest" cell in the boundary and subtract it from
				// the current cell.
				sym = Set();
				youngest = Boundary[nextColumnAdded[Mod2youngestOf(cell)]];
				std::set_symmetric_difference(cell.begin(), cell.end(), youngest.begin(), youngest.end(), std::inserter(sym, sym.begin()));
				cell = sym;
			}
			// Check whether we've eliminated the column. For some god damn reason
			// we have to re-set the entry of the Boundary??? Why?????? Scope??? wtf
			if (!cell.empty()) {
				Boundary[j] = cell;
				nextColumnAdded[Mod2youngestOf(cell)] = j;
				Boundary[Mod2youngestOf(cell)] = Set();
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
		if (nextColumnAdded[*it] == 0) essential.insert(*it);
	}

	return essential;
}


Set ComputePercolationEvents(
		Vector boundary, Vector filtration, int homology, int p, Vector breaks
	) {
	// Disable the Z/2Z check for now.	
	// if (p == 2) { return Mod2ComputePercolationEvents(boundary, filtration, homology, breaks); }

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

	Field F(p);
	Vector _boundary = ReindexBoundaryMatrix(boundary, filtration, homology, breaks);
	BoundaryMatrix Boundary = FillBoundaryMatrix(_boundary, F, cellCount);

	Vector nextColumnAdded = Vector(cellCount, 0);
	Column cell, youngest;
	Set erase, marked;
	int face, high, numBreaks = breaks.size();

	Field::Element q, r, s, inv, prod, result;
	F.init(inv);
	F.init(prod);
	F.init(result);

	marked = Set();

	for (int d = topDimension; d > homology-1; d--) {

		high = (d+1 >= numBreaks ? cellCount : breaks[d+1]);

		for (int j = breaks[d]; j < high; j++) {
			// If we're of the wrong dimension, keep going.
			// if (dim(j, breaks) != d) { continue; }
			cell = Boundary[j];

			while (!cell.empty() && nextColumnAdded[youngestOf(cell)] != 0) {
				// Keep track of keys to erase.
				erase = Set();

				// Get the "youngest" cell in the boundary and subtract it from
				// the current cell.
				youngest = Boundary[nextColumnAdded[youngestOf(cell)]];

				// Get the multiplicative inverse of the coefficient and do
				// arithmetic over the row.
				q = youngest[youngestOf(cell)];
				F.inv(inv, q);

				for (auto it=youngest.begin(); it != youngest.end(); ++it) {
					face = it->first;
					s = it->second;

					// Take the product of inv(q) with the entry of this row;
					// if this is row youngestOf(cell), then this product is 1
					// (it's a pivot). If the current cell shares this face, do
					// the subtraction; otherwise, just add a new coefficient.
					F.mul(prod, inv, s);

					if (cell.count(face) > 0) {
						F.sub(result, cell[face], prod);

						if (F.isZero(result)) {
							cell.erase(face);
						} else {
							cell[face] = result;
						}
					} else {
						F.neg(result, prod);
						cell[face] = result;
					}
				}
			}
			// Check whether we've eliminated the column. For some god damn reason
			// we have to re-set the entry of the Boundary??? Why?????? Scope??? wtf
			if (!cell.empty()) {
				Boundary[j] = cell;
				nextColumnAdded[youngestOf(cell)] = j;
				Boundary[youngestOf(cell)] = Column();
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
		if (nextColumnAdded[*it] == 0) essential.insert(*it);
	}

	return essential;
}
