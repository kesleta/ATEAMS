
#include <fstream>
#include <iostream>
#include <string>
#include <ATEAMS/common.h>
#include <ATEAMS/LinBoxMethods.h>

using namespace std;


Index matrix(string filename) {
	Index matrix;
	ifstream file(filename);
	string line;

	while(getline(file, line)) {
		matrix.push_back((INDEXTYPE)stoi(line));
	}

	return matrix;
}

BoundaryMatrix FillBoundaryMatrix(Index boundary, int cellCount) {
	BoundaryMatrix Boundary(cellCount, Column());
	int row, column, q;

	for (int t = 0; t < boundary.size(); t+=3) {
		row = boundary[t];
		column = boundary[t+1];
		q = boundary[t+2];

		if (q != 0) { Boundary[column][row] = (DATATYPE)q; }
	}

	return Boundary;
}

int main() {
	Index boundary = matrix("bd.txt");
	Index filtration = matrix("filtration.txt");
	Index breaks = {0, 10000, 50000, 110000, 150000};

	int field = 3;
	int homology = 2;
	int topDimension = 4;
	int cellCount = 160000;

	BoundaryMatrix Boundary = FillBoundaryMatrix(boundary, cellCount);
	set<int> essential = ZpComputePercolationEvents(field, Boundary, breaks, cellCount);

	for (int t=0; t < 100; t++) {
		set<int> essential = ZpComputePercolationEvents(field, Boundary, breaks, cellCount);
		
		// cout << "{ ";
		// for (auto it=essential.begin(); it != essential.end(); it++) {
		// 	cout << *it << " ";
		// }
		// cout << "}" << endl;
		// cout << endl;
	}

	return 0;
}

