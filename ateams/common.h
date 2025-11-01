

#ifndef ATEAMS_COMMON_H
#define ATEAMS_COMMON_H

#include <set>
#include <vector>
#include <map>

#include <unistd.h>
#include <stdio.h>

using namespace std;

// Type for indices and type for data; for now, we use unsigned 32-bit ints for
// indices (since we're never negative) and signed chars for data (since our
// fields are never that big).
typedef int32_t INDEXTYPE;
typedef char DATATYPE;

// Standard sets, indices
typedef set<INDEXTYPE> Set;
typedef map<INDEXTYPE,INDEXTYPE> Map;
typedef vector<INDEXTYPE> Index;

// Arithmetic.
typedef vector<DATATYPE> Lookup;
typedef vector<Lookup> Table;

// Persistence pairs; these are literally just pairs of ints.
typedef vector<Index> PersistencePairs;

// Matrices.
typedef map<INDEXTYPE,DATATYPE> Column;
typedef vector<Column> BoundaryMatrix;
typedef vector<Index> FlatBoundaryMatrix;

// Utility functions for suppressing output to stderr (specifically because SpaSM)
// outputs diagnostic information at every step!).
int _suppress() {
	fflush(stderr);
	int fd = dup(STDERR_FILENO);
	freopen("/dev/null", "w", stderr);
	return fd;
}

void _resume(int fd) {
	fflush(stderr);
	dup2(fd, fileno(stderr));
	close(fd);
}

#endif

