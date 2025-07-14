
#include <set>
#include <vector>
using namespace std;

typedef set<int> Set;
typedef vector<int> Index;
typedef vector<char> Lookup;
typedef vector<Lookup> Table;
typedef map<int,char> Column;
typedef vector<Column> BoundaryMatrix;

Index LanczosKernelSample(Index coboundary, int M, int N, int p, int maxTries);
Set ComputePercolationEvents(Table addition, Table multiplication, Lookup negation, Lookup inversion, BoundaryMatrix Boundary, Index breaks, int cellCount, int topDimension, int homology);
Set ZpComputePercolationEvents(int field, BoundaryMatrix Boundary, Index breaks, int cellCount, int topDimension, int homology);
