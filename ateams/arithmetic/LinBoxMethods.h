

#include <ATEAMS/common.h>

Index ReducedKernelSample(Index coboundary, int M, int N, int p);
Set ComputePercolationEvents(Table addition, Table multiplication, Lookup negation, Lookup inversion, BoundaryMatrix Boundary, Index breaks, int cellCount);
Set LinearComputePercolationEvents(int field, Lookup addition, Lookup multiplication, Lookup negation, Lookup inversion, BoundaryMatrix Boundary, Index breaks, int cellCount, int dimension);
Set ZpComputePercolationEvents(int field, BoundaryMatrix Boundary, Index breaks, int cellCount);
