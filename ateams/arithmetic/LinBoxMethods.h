
using namespace std;

std::vector<int> LanczosKernelSample(std::vector<int> coboundary, int M, int N, int p, int maxTries);
std::set<int> ComputePercolationEvents(std::vector<int> boundary, std::vector<int> filtration, int homology, int p, vector<int> breaks);
