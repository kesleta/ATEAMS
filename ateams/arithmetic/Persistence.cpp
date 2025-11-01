
#include <phat/compute_persistence_pairs.h>
#include <phat/boundary_matrix.h>
#include <phat/representations/default_representations.h>
#include <phat/algorithms/twist_reduction.h>
#include <phat/algorithms/standard_reduction.h>
#include <phat/algorithms/row_reduction.h>
#include <phat/algorithms/chunk_reduction.h>
#include <phat/algorithms/spectral_sequence_reduction.h>
#include <phat/algorithms/swap_twist_reduction.h>
#include <phat/algorithms/exhaustive_compress_reduction.h>
#include <phat/algorithms/lazy_retrospective_reduction.h>
#include <phat/helpers/dualize.h>

#include "Persistence.h"

using namespace std;

/*
#################################
##### PERSISTENCE WITH PHAT #####
#################################
*/
typedef vector<phat::index> PHATColumn;
typedef phat::boundary_matrix<phat::bit_tree_pivot_column> PHATBoundaryMatrix;
typedef phat::persistence_pairs Pairs;
typedef phat::twist_reduction Twist;


Index ReindexBoundaryMatrix(Index &boundary, Index filtration, int homology, Index breaks) {
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


FlatBoundaryMatrix Flatten(Index full, int columns, Index breaks) {
	FlatBoundaryMatrix boundary = FlatBoundaryMatrix(columns, Index());
	int row, column;

	for (int t=0; t < full.size(); t+=3) {
		row = full[t];
		column = full[t+1];

		// If the column corresponds to a cell that isn't a vertex, we insert the
		// row into the column. Afterwards, we'll have to go through and sort
		// each of the columns as well.
		if (breaks[1] <= column) { boundary[column].push_back(row); }
	}

	// Sort column entries (presumably so this makes building the bitpacks easier
	// later?).
	for (int i=0; i < columns; i++) {
		sort(boundary[i].begin(), boundary[i].end());
	}

	return boundary;
}


FlatBoundaryMatrix ReindexAndFlatten(Index _boundary, Index filtration, int homology, Index breaks) {
	_boundary = ReindexBoundaryMatrix(_boundary, filtration, homology, breaks);
	return Flatten(_boundary, filtration.size(), breaks);
}


PersistencePairs PHATComputePersistencePairs(Index _boundary, Index filtration, int homology, Index breaks) {
	// Build out the boundary matrix.
	Pairs pairs;
	FlatBoundaryMatrix boundary = ReindexAndFlatten(_boundary, filtration, homology, breaks);
	PHATBoundaryMatrix Boundary;
	PHATColumn faces;
	int numFaces, numPHATColumns = boundary.size();

	// Count the number of columns.
	Boundary.set_num_cols(numPHATColumns);

	for (int t=0; t < numPHATColumns; t++) {
		// Fill in the faces for this column. Should be fine dealing with empty
		// columns (i.e. vertices).
		numFaces = boundary[t].size();

		for (int j=0; j < numFaces; j++) { faces.push_back(boundary[t][j]); }

		Boundary.set_dim(t, numFaces*2);
		Boundary.set_col(t, faces);
		faces.clear();
	}

	// Compute the persistence pairs and populate a Vector to return to the
	// user.
	phat::compute_persistence_pairs<Twist>(pairs, Boundary);
	pairs.sort();

	PersistencePairs Pairs(pairs.get_num_pairs(), Index(2));

	for (phat::index i=0; i < pairs.get_num_pairs(); i++) {
		Pairs[i][0] = pairs.get_pair(i).first;
		Pairs[i][1] = pairs.get_pair(i).second;
	}

	return Pairs;
}


/*
#########################################
##### PERSISTENCE WITH LINBOX/SPASM #####
#########################################
*/
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/ring/modular.h>

typedef Givaro::Modular<int> Zp;
typedef LinBox::SparseMatrix<Zp, LinBox::SparseMatrixFormat::SparseSeq> ZpMatrix;
typedef LinBox::DenseVector<Zp> ZpVector;


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

