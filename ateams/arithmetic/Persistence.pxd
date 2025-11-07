
from ..common cimport Index, PersistencePairs, Table, Lookup, BoundaryMatrix, Set, Column, Bases, Basis, bool

cdef extern from "Persistence.h":
	# PHAT solutions.
	PersistencePairs PHATComputePersistencePairs(Index boundary, Index filtration, int homology, Index breaks) noexcept

	# Homebrew solutions.
	Set ComputePercolationEvents(Table addition, Table multiplication, Lookup negation, Lookup inversion, BoundaryMatrix Boundary, Index breaks, int cellCount) noexcept
	Set LinearComputePercolationEvents(int field, Lookup addition, Lookup multiplication, Lookup negation, Lookup inversion, BoundaryMatrix Boundary, Index breaks, int cellCount, int dimension) noexcept
	Set ZpComputePercolationEvents(int field, BoundaryMatrix Boundary, Index breaks, int cellCount) noexcept

	# SpaSM solutions.
	Bases LinearComputeBases(int field, Lookup addition, Lookup multiplication, Lookup negation, Lookup inversion, BoundaryMatrix Boundary, Index breaks, int cellCount, int dimension);
	Basis ComputeCobasis(Basis combined, int M, int N, int rank, int characteristic, bool verbose);
	Set RankComputePercolationEvents(BoundaryMatrix augmented, int M, int N, int basisrank, int p, bool verbose);
	Set RankComputePercolationEvents(BoundaryMatrix augmented, int M, int N, int basisrank, int p, int stop, bool verbose);

	# Deprecated.
	Set SolveComputePercolationEvents(BoundaryMatrix coboundary, Basis cobasis, int M, int N, int basisrank, int p);

	# SparseRREF solutions.
	Set SRankComputePercolationEvents(BoundaryMatrix augmented, int M, int N, int basisrank, int p, bool verbose);

