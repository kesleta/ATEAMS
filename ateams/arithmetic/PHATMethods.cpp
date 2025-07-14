
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

#include "PHATMethods.h"

using namespace std;


typedef vector<phat::index> Column;
typedef phat::boundary_matrix<phat::bit_tree_pivot_column> BoundaryMatrix;
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
	BoundaryMatrix Boundary;
	Column faces;
	int numFaces, numColumns = boundary.size();

	// Count the number of columns.
	Boundary.set_num_cols(numColumns);

	for (int t=0; t < numColumns; t++) {
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
