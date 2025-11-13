
#ifndef ATEAMS_PERSISTENCE_H
#define ATEAMS_PERSISTENCE_H

#include <ATEAMS/common.h>

// PHAT solutions.
PersistencePairs PHATComputePersistencePairs(Index boundary, Index filtration, int homology, Index breaks);

// Homebrew solutions.
Set ComputePercolationEvents(Table addition, Table multiplication, Lookup negation, Lookup inversion, BoundaryMatrix Boundary, Index breaks, int cellCount);
Set LinearComputePercolationEvents(int field, Lookup addition, Lookup multiplication, Lookup negation, Lookup inversion, BoundaryMatrix Boundary, Index breaks, int cellCount, int dimension);
Set ZpComputePercolationEvents(int field, BoundaryMatrix Boundary, Index breaks, int cellCount);

// SpaSM solutions.
Bases LinearComputeBases(int field, Lookup addition, Lookup multiplication, Lookup negation, Lookup inversion, BoundaryMatrix Boundary, Index breaks, int cellCount, int dimension);
Basis ComputeCobasis(Basis combined, int M, int N, int rank, int characteristic, bool verbose=true);
Set RankComputePercolationEvents(BoundaryMatrix augmented, int M, int N, int basisrank, int p, bool verbose=true);
Set RankComputePercolationEvents(BoundaryMatrix augmented, int M, int N, int basisrank, int p, int stop, bool verbose=true);

// Deprecated.
Set SolveComputePercolationEvents(BoundaryMatrix coboundary, Basis cobasis, int M, int N, int basisrank, int p);

// SparseRREF solutions.
Set SRankComputePercolationEvents(BoundaryMatrix augmented, int M, int N, int basisrank, int p, bool verbose=true);
Set SRankComputePercolationEvents(BoundaryMatrix augmented, int M, int N, int basisrank, int p, int stop, bool verbose=true);

#endif
