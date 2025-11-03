
#include <SparseRREF/sparse_mat.h>

#include "Sampling.h"
#include <iostream>

using namespace std;

typedef SparseRREF::sparse_mat<data_t, index_t> ZpMatrix;
typedef SparseRREF::sparse_vec<data_t, index_t> ZpVector;
typedef SparseRREF::field_t Zp;


Index randomVector(int N, int p) {
	Index r(N, 0);
	for (int i=0; i<N; i++) r[i] = (data_t)(random()%p);
	return r;
}

template <typename Matrix, typename Vector, typename Field>
Index RandomLinearCombination(Matrix K, Field GFp, int p) {
	Matrix KT = K.transpose();
	Index coeff = randomVector(KT.nrow, p);
	Index rlc = Index(KT.ncol, 0);

	cerr << "[Sampling] building random kernel element...";
	for (int i=0; i<KT.nrow; i++) {
		for (auto [j,q] : KT.rows[i]) {
			rlc[j] = SparseRREF::scalar_add(rlc[j], SparseRREF::scalar_mul(coeff[i], q, GFp), GFp);
		}
	}
	cerr << " done." << endl;

	return rlc;
}


template <typename Matrix>
Matrix SparseMatrixFill(Index coboundary, int M, int N, int p) {
	Matrix A(M, N);
	int row, col, _val, val;

	cerr << "[Sampling] building sparse coboundary matrix...";
	for (int t=0; t<coboundary.size(); t+=3) {
		// Pull out the values from the coboundary matrix.
		row = coboundary[t];
		col = coboundary[t+1];
		_val = coboundary[t+2];
		val = _val < 0 ? (p-1) : 1;

		// Insert them into the matrix.
		ZpVector& matrow = A.rows[row];
		matrow.push_back((index_t)col, (data_t)val);
	}

	A.compress();
	cerr << " done." << endl;
	return A;
}


Index ZpReducedKernelSample(Index coboundary, int M, int N, int p) {
	// Construct field, matrix, and find a basis for the kernel.
	const Zp GFp(SparseRREF::FIELD_Fp, p);
	ZpMatrix A = SparseMatrixFill<ZpMatrix>(coboundary, M, N, p);

	SparseRREF::rref_option_t opt;
	opt->pool.reset();
	opt->method = 0;

	std::vector<std::vector<SparseRREF::pivot_t<index_t>>> pivots;
	cerr << "[Sampling] RREFing...";
	pivots = SparseRREF::sparse_mat_rref<data_t, index_t>(A, GFp, opt);
	cerr << " done." << endl;

	cerr << "[Sampling] computing kernel (RREF)...";
	auto K = SparseRREF::sparse_mat_rref_kernel<data_t, index_t>(A, pivots, GFp, opt);
	cerr << " done." << endl;

	opt->abort = true;
	Flint::clear_cache();
	return RandomLinearCombination<ZpMatrix, ZpVector>(K, GFp, p);
}

Index ReducedKernelSample(Index coboundary, int M, int N, int p) {
	cerr << endl;
	return ZpReducedKernelSample(coboundary, M, N, p);
}
