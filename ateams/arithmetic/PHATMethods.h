
#include <vector>
#include <map>

using namespace std;

typedef vector<int> Index;
typedef vector<Index> PersistencePairs;
typedef vector<Index> FlatBoundaryMatrix;
typedef map<int,int> Map;


PersistencePairs PHATComputePersistencePairs(Index boundary, Index filtration, int homology, Index breaks);
