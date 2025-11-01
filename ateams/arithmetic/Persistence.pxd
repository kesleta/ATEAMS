
from ..common cimport Index, PersistencePairs, Table, Lookup, BoundaryMatrix, Set

cdef extern from "Persistence.h":
	PersistencePairs PHATComputePersistencePairs(Index boundary, Index filtration, int homology, Index breaks) noexcept
	Set ComputePercolationEvents(Table addition, Table multiplication, Lookup negation, Lookup inversion, BoundaryMatrix Boundary, Index breaks, int cellCount) noexcept
	Set LinearComputePercolationEvents(int field, Lookup addition, Lookup multiplication, Lookup negation, Lookup inversion, BoundaryMatrix Boundary, Index breaks, int cellCount, int dimension) noexcept
	Set ZpComputePercolationEvents(int field, BoundaryMatrix Boundary, Index breaks, int cellCount) noexcept

