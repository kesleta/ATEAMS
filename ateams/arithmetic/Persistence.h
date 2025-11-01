
#ifndef ATEAMS_PERSISTENCE_H
#define ATEAMS_PERSISTENCE_H

#include <ATEAMS/common.h>

PersistencePairs PHATComputePersistencePairs(Index boundary, Index filtration, int homology, Index breaks);
Set ComputePercolationEvents(Table addition, Table multiplication, Lookup negation, Lookup inversion, BoundaryMatrix Boundary, Index breaks, int cellCount);
Set LinearComputePercolationEvents(int field, Lookup addition, Lookup multiplication, Lookup negation, Lookup inversion, BoundaryMatrix Boundary, Index breaks, int cellCount, int dimension);
Set ZpComputePercolationEvents(int field, BoundaryMatrix Boundary, Index breaks, int cellCount);

#endif
