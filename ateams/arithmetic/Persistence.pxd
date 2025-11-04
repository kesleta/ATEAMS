
from ..common cimport Index, PersistencePairs, Table, Lookup, BoundaryMatrix, Set, Column, Bases, Basis

cdef extern from "Persistence.h":
	PersistencePairs PHATComputePersistencePairs(Index boundary, Index filtration, int homology, Index breaks) noexcept
	Set ComputePercolationEvents(Table addition, Table multiplication, Lookup negation, Lookup inversion, BoundaryMatrix Boundary, Index breaks, int cellCount) noexcept
	Set LinearComputePercolationEvents(int field, Lookup addition, Lookup multiplication, Lookup negation, Lookup inversion, BoundaryMatrix Boundary, Index breaks, int cellCount, int dimension) noexcept
	Set ZpComputePercolationEvents(int field, BoundaryMatrix Boundary, Index breaks, int cellCount) noexcept

	Bases LinearComputeBases(int field, Lookup addition, Lookup multiplication, Lookup negation, Lookup inversion, BoundaryMatrix Boundary, Index breaks, int cellCount, int dimension);
	Basis ComputeCobasis(Basis combined, int M, int N, int rank, int characteristic);
	Set RankComputePercolationEvents(BoundaryMatrix augmented, int M, int N, int basisrank, int p);
	Set SolveComputePercolationEvents(BoundaryMatrix coboundary, Basis cobasis, int M, int N, int basisrank, int p);

