
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/ring/modular.h>

extern "C" {
	#include <spasm/spasm.h>
}

#include "LinBoxMethods.h"
#include <iostream>

using namespace std;

typedef Givaro::Modular<int> Zp;
typedef LinBox::SparseMatrix<Zp, LinBox::SparseMatrixFormat::SparseSeq> ZpMatrix;
typedef LinBox::DenseVector<Zp> ZpVector;


template <typename Vector, typename Field>
Vector randomVector(int N, Field F) {
	Vector r(F, N);

	for (int i=0; i<N; i++) {
		// Here we're presuming the field has prime --- not prime *power* --- order.
		typename Field::Element q;
		F.init(q, random()%(int)F.characteristic());
		r.setEntry(i, q);
	}

	return r;
}


template <typename Matrix, typename Field>
Matrix MatrixFromSpaSM(const struct spasm_csr *A, Field F) {
	// Get the dimensions of the spasm_csr matrix. NOTE: SpaSM uses n and m to
	// refer to the number of rows and columns of a matrix, contrary to convention.
	// What's below is not a typo!
	int M = A->n;
	int N = A->m;
	Matrix B(F, M, N);

	// Get pointers to column indices and rows.
	const int *columnIndices = A->j;
	const i64 *rows = A->p;
	const spasm_ZZp *values = A->x;

	// Iterate through the CSR, adding entries to B.
	for (int row=0; row<M; row++) {
		for (i64 px=rows[row]; px<rows[row+1]; px++) {
			i64 x = values[px];

			typename Field::Element q;
			F.init(q, x);
			B.setEntry(row, columnIndices[px], q);
		}
	}

	B.finalize();
	return B;
}


template <typename Vector>
Index populate(Vector X) {
	// Populate a vector with entries from a LinBox vector.
	Index x(X.size());
	for (size_t k = 0; k < X.size(); k++) { x[k] = X.getEntry(k); }

	return x;
}


Index ReducedKernelSample(Index coboundary, int M, int N, int p) {
	// Construct the finite field and construct the matrix. If we're over Z/2Z,
	// use specialized matrices for our operations.
	typedef ZpMatrix Matrix;
	typedef ZpVector Vector;
	typedef Zp Field;

	// Construct the field and vectors; allocate a SpaSM matrix.
	Field F(p);
	ZpVector X(F, N);

	int _sup = _suppress();
	struct spasm_triplet *T = spasm_triplet_alloc(M, N, 1, p, p!=-1);
	for (int t=0; t<coboundary.size(); t+=3) { spasm_add_entry(T, coboundary[t], coboundary[t+1], coboundary[t+2]); }

	struct spasm_csr *A = spasm_compress(T);
	spasm_triplet_free(T);

	// Echelonize A and compute a basis for the kernel.
	struct spasm_lu *fact = spasm_echelonize(A, NULL);
	spasm_csr_free(A);

	const struct spasm_csr *K = spasm_kernel(fact);
	Matrix basis = MatrixFromSpaSM<Matrix, Field>(K, F);
	_resume(_sup);
	
	// Take the kernel basis, read it into a LinBox matrix, generate some random
	// coefficients, then take a linear combination of the rows of the basis to
	// get a random vector in the kernel.
	Vector combination = randomVector<Vector, Field>(basis.rowdim(), F);
	Vector assignment(F, N);
	assignment = basis.applyTranspose(assignment, combination);

	return populate(assignment);
}


template <typename BalancedStorage>
int youngestOf(BalancedStorage column) {
	// Gets the "youngest" (largest-indexed) cell in the column.
	return column.rbegin()->first;
}


typedef map<int,Zp::Element> ZpColumn;
typedef vector<ZpColumn> ZpBoundaryMatrix;


ZpBoundaryMatrix ZpFillBoundaryMatrix(BoundaryMatrix Boundary, Zp Field) {
	ZpBoundaryMatrix B(Boundary.size());
	Zp::Element q;
	Column column;
	ZpColumn pcolumn;
	int face;

	for (int t = 0; t < Boundary.size(); t++) {
		column = Boundary[t];
		pcolumn = ZpColumn();

		for (auto it=column.begin(); it != column.end(); it++) {
			face = it->first;
			Field.init(q, it->second);
			pcolumn[face] = q;
		}

		B[t] = pcolumn;
	}

	return B;
}


Set ZpComputePercolationEvents(int field, BoundaryMatrix _boundary, Index breaks, int cellCount) {
	Zp F(field);
	ZpBoundaryMatrix Boundary = ZpFillBoundaryMatrix(_boundary, F);

	Index nextColumnAdded = Index(cellCount, 0);
	ZpColumn cell, youngest;
	Set marked = Set();
	int face, high, numBreaks = breaks.size();

	Zp::Element q, r, s, inv, neg, prod, result;
	F.init(q);
	F.init(r);
	F.init(s);
	// char q, r, s, inv, prod, result;

	for (int d = numBreaks-1; d > 0; d--) {

		high = (d+1 >= numBreaks ? cellCount : breaks[d+1]);

		for (int j = breaks[d]; j < high; j++) {
			// If we're of the wrong dimension, keep going.
			// if (dim(j, breaks) != d) { continue; }
			cell = Boundary[j];

			while (!cell.empty() && nextColumnAdded[youngestOf(cell)] != 0) {
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
					// prod = multiplication[inv][s];

					if (cell.count(face) > 0) {
						F.sub(result, cell[face], prod);
						// result = addition[cell[face]][negation[prod]];

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
				Boundary[youngestOf(cell)] = ZpColumn();
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
		Table addition, Table multiplication, Lookup negation, Lookup inversion,
		BoundaryMatrix Boundary, Index breaks, int cellCount
	) {
	Index nextColumnAdded = Index(cellCount, 0);
	Column cell, youngest;
	Set marked = Set();
	int face, high, numBreaks = breaks.size();

	char q, r, s, inv, prod, result;

	for (int d = numBreaks-1; d > 0; d--) {

		high = (d+1 >= numBreaks ? cellCount : breaks[d+1]);

		for (int j = breaks[d]; j < high; j++) {
			// If we're of the wrong dimension, keep going.
			// if (dim(j, breaks) != d) { continue; }
			cell = Boundary[j];

			while (!cell.empty() && nextColumnAdded[youngestOf(cell)] != 0) {
				// Get the "youngest" cell in the boundary and subtract it from
				// the current cell.
				youngest = Boundary[nextColumnAdded[youngestOf(cell)]];

				// Get the multiplicative inverse of the coefficient and do
				// arithmetic over the row.
				q = youngest[youngestOf(cell)];
				inv = inversion[q];

				for (auto it=youngest.begin(); it != youngest.end(); ++it) {
					face = it->first;
					s = it->second;

					// Take the product of inv(q) with the entry of this row;
					// if this is row youngestOf(cell), then this product is 1
					// (it's a pivot). If the current cell shares this face, do
					// the subtraction; otherwise, just add a new coefficient.
					prod = multiplication[inv][s];

					if (cell.count(face) > 0) {
						result = addition[cell[face]][negation[prod]];

						if (result < 1) {
							cell.erase(face);
						} else {
							cell[face] = result;
						}
					} else {
						cell[face] = negation[prod];
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

typedef function<DATATYPE(DATATYPE, DATATYPE)> binop;
typedef function<DATATYPE(DATATYPE)> unop;


binop _add(Lookup addition, DATATYPE p) {
	return [addition, p](DATATYPE a, DATATYPE b) { return addition[a*p + b]; };
}

binop _mult(Lookup multiplication, DATATYPE p) {
	return [multiplication, p](DATATYPE a, DATATYPE b) { return multiplication[a*p + b]; };
}

unop _neg(Lookup negation) {
	return [negation](DATATYPE a) { return negation[a]; };
}

unop _inv(Lookup inversion) {
	return [inversion](DATATYPE a) { return inversion[a]; };
}


Set LinearComputePercolationEvents(
		int field, Lookup addition, Lookup multiplication, Lookup negation, Lookup inversion,
		BoundaryMatrix Boundary, Index breaks, int cellCount, int dimension
	) {
	Index nextColumnAdded = Index(cellCount, 0);
	Column cell, youngest;
	Set marked = Set();
	int face, high, numBreaks = breaks.size();

	binop add = _add(addition, field);
	binop multiply = _mult(multiplication, field);
	unop negate = _neg(negation);
	unop invert = _inv(inversion);

	char q, r, s, inv, prod, result;

	for (int d = numBreaks-1; d > dimension-1; d--) {

		high = (d+1 >= numBreaks ? cellCount : breaks[d+1]);

		for (int j = breaks[d]; j < high; j++) {
			// If we're of the wrong dimension, keep going.
			// if (dim(j, breaks) != d) { continue; }
			cell = Boundary[j];

			while (!cell.empty() && nextColumnAdded[youngestOf(cell)] != 0) {
				// Get the "youngest" cell in the boundary and subtract it from
				// the current cell.
				youngest = Boundary[nextColumnAdded[youngestOf(cell)]];

				// Get the multiplicative inverse of the coefficient and do
				// arithmetic over the row.
				q = youngest[youngestOf(cell)];
				inv = invert(q);

				for (auto it=youngest.begin(); it != youngest.end(); ++it) {
					face = it->first;
					s = it->second;

					// Take the product of inv(q) with the entry of this row;
					// if this is row youngestOf(cell), then this product is 1
					// (it's a pivot). If the current cell shares this face, do
					// the subtraction; otherwise, just add a new coefficient.
					prod = multiply(inv, s);
					// prod = multiplication[inv][s];

					if (cell.count(face) > 0) {
						result = add(cell[face], negate(prod));
						// result = addition[cell[face]][negation[prod]];

						if (result < 1) {
							cell.erase(face);
						} else {
							cell[face] = result;
						}
					} else {
						cell[face] = negate(prod);
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
