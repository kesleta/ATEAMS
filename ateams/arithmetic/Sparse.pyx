
# cython: language_level=3str, initializedcheck=False, c_api_binop_methods=True, nonecheck=False, profile=True, cdivision=True
# cython: binding=True, linetrace=True
# cython: boundscheck=False, wraparound=False
# define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
# distutils: language=c++

import numpy as np
cimport numpy as np

from .common cimport FFINT, FLAT, TABLE

from cython.parallel cimport prange
from openmp cimport omp_get_thread_num as Thread
from libcpp cimport bool
from libcpp.unordered_set cimport unordered_set as Set
from libcpp.vector cimport vector as Vector
from libcpp.unordered_map cimport unordered_map as Map


cdef Set[int] _intersect = Set[int]();
cdef Set[int] _union = Set[int]();


cdef Set[int] Difference(Set[int] P, Set[int] Q) noexcept:
	cdef Set[int] D = Set[int](P);
	cdef int s;

	for s in P:
		if Q.contains(s): D.erase(s);

	return D


cdef Set[int] Union(Set[int] P, Set[int] Q) noexcept:
	cdef Set[int] S = Set[int](P);
	cdef int s;

	for s in Q: S.insert(s);

	return S;


cdef Set[int] Intersection(Set[int] P, Set[int] Q) noexcept:
	cdef Set[int] C, T, S;
	cdef int s;

	_intersect.clear();

	# Only iterate over the smaller set.
	if P.size() < Q.size():
		T = P;
		S = Q;
	else:
		T = Q;
		S = P;

	for s in T:
		if S.contains(s): _intersect.insert(s)

	return _intersect


cdef class Matrix:
	def __cinit__(
			self,
			TABLE A,
			TABLE addition,
			FLAT negation,
			TABLE multiplication,
			FLAT inverses,
			bool parallel,
			int minBlockSize,
			int maxBlockSize,
			int cores
		):
		# Initialize addition and multiplication tables.
		self.addition = addition;
		self.negation = negation;
		self.multiplication = multiplication;
		self.inverses = inverses;
		self.data = A;

		# Initialize the `.shape` property.
		self.shape = Vector[int]();
		self.shape.push_back(A.shape[0]);
		self.shape.push_back(A.shape[1]);

		# If we want synchronous matrix operations, set the number of cores,
		# the minimum block size, the threading schedule, and the block schema.
		self.parallel = parallel;
		self.cores = cores;
		self.minBlockSize = minBlockSize;
		self.maxBlockSize = maxBlockSize;
		self.blockSchema = Vector[Vector[int]]();

		# Initialize the `.rows` and `.columns` properties.
		self.columns = Map[int, Set[int]]();
		self._initializeColumns();

		# Keep track of the lowest zero row.
		self.zero = -1


	cdef Vector[Vector[int]] recomputeBlockSchema(self, int start, int stop) noexcept:
		cdef int block, blocks, MINCOL, MAXCOL, b, _q;
		cdef Vector[int] BLOCK;

		self.blockSchema.clear();
		
		# Compute the block size *of the current submatrix*, then schedule vertical
		# block computation appropriately: we have bounds on the "block size,"
		# so no thread is given too much (or too little) work to do. Choosing the
		# block sizes may be more of an art than a science, though.
		_q = (stop-start)//self.cores;
		block = min(max(self.minBlockSize, _q), min(self.maxBlockSize, _q))
		blocks = ((stop-start)//block)+1
		BLOCKS = Vector[Vector[int]]();

		for b in range(blocks):
			MINCOL = start+b*block;
			MAXCOL = start+(b+1)*block if start+(b+1)*block <= stop else stop
			
			if MINCOL >= stop: break

			BLOCK = Vector[int]();
			BLOCK.push_back(MINCOL);
			BLOCK.push_back(MAXCOL);
			self.blockSchema.push_back(BLOCK);

		return self.blockSchema;
	

	cdef void _initializeColumns(self) noexcept:
		"""
		Initializes the `.rows` and `.columns` data structures.
		"""
		cdef int i, j, M, N;
		cdef Set[int] columns;

		M = self.shape[0];
		N = self.shape[1];

		for i in range(M):
			columns = Set[int]();

			for j in range(N):
				if self.data[i,j] > 0:
					columns.insert(j);
			
			self.columns[i] = columns;
	
	
	cdef void SwapRows(self, int i, int j):
		"""
		Swaps rows `i` and `j`.
		"""
		# Swap entries in the database, re-scanning rows first.
		cdef int b, c;
		cdef FFINT t;
		
		if self.parallel:
			self.ScanRow(i);
			self.ScanRow(j);
		
		cdef Set[int] COLUMNS = Union(self.columns[i], self.columns[j]);

		for c in COLUMNS:
			t = self.data[i,c];
			self.data[i,c] = self.data[j,c];
			self.data[j,c] = t;

		# Swap the columns (just by swapping keys).
		cdef Set[int] colswap = self.columns[i];
		self.columns[i] = self.columns[j];
		self.columns[j] = colswap;


	cdef void ParallelAddRows(self, int i, int j, int start, int stop, FFINT ratio) noexcept nogil:
		"""
		Adds rows `i` and `j`.
		"""
		# First, perform all the addition, then scan over both rows to delete fill
		# (extra added zeros). WARNING: this operation is *not commutative*! If
		# we're adding row `i` to row `j`, then we only add *nonzero elements*
		# of row `i` to row `j`.
		cdef int c, t;
		cdef FFINT p, q, mult;

		for c in self.columns[i]:
			if c < start or c >= stop: continue

			p = self.data[i,c]
			mult = self.multiplication[p, ratio]
			self.data[j,c] = self.addition[mult,self.data[j,c]]
	

	cdef void AddRows(self, int i, int j, int start, int stop, FFINT ratio) noexcept:
		"""
		Adds rows `i` and `j`.
		"""
		# First, perform all the addition, then scan over both rows to delete fill
		# (extra added zeros). WARNING: this operation is *not commutative*! If
		# we're adding row `i` to row `j`, then we only add *nonzero elements*
		# of row `i` to row `j`.
		cdef int c;
		cdef FFINT p, q, mult;

		for c in self.columns[i]:
			if c < start or c >= stop: continue

			p = self.data[i,c]
			mult = self.multiplication[p, ratio]
			self.data[j,c] = self.addition[mult,self.data[j,c]]

			# Update columns.
			if self.data[j,c] < 1: self.columns[j].erase(c)
			else: self.columns[j].insert(c)


	cdef void MultiplyRow(self, int i, FFINT q) noexcept:
		"""
		Multiplies row `i` by `q`.
		"""
		cdef int c;
		cdef FFINT p;

		for c in self.columns[i]:
			p = self.data[i,c];
			self.data[i,c] = self.multiplication[q,p];

	
	cdef void ScanRow(self, int i):
		cdef int c;
		
		for c in range(self.shape[1]):
			if self.data[i,c] < 1: self.columns[i].erase(c);
			else: self.columns[i].insert(c);
		

	cdef int PivotRow(self, int c, int pivots) noexcept:
		"""
		Report the first nonzero entry in column `c`; if no such entry exists,
		return -1.
		"""
		cdef int i, p;
		p = -1;

		for i in range(self.shape[0]):
			if self.parallel:
				if self.data[i,c] > 0:
					self.columns[i].insert(c);
					p = i if p < 0 and i > pivots-1 else p
				else:
					self.columns[i].erase(c);
			else:
				if self.data[i,c] > 0 and i > pivots-1:
					p = i;
					break;

		return p

	
	cdef int HighestZeroRow(self, int AUGMENT=-1) noexcept:
		"""
		Report the first lowest-indexed zero row; if no such row exists,
		report -1. If nonnegative AUGMENT is passed, we check for the first row
		with AUGMENT leading zeros.
		"""
		# Catch for optional arguments.
		if AUGMENT < 0: AUGMENT = self.shape[1]
		return self.zero if self.zero <= AUGMENT else -1
	

	cdef TABLE ToArray(self) noexcept: return self.data
	

	cdef void RREF(self, int AUGMENT=-1) noexcept:
		"""
		Computes the RREF of this matrix.
		"""
		# Catch for optional arguments.
		if AUGMENT < 0: AUGMENT = self.shape[1]

		cdef Vector[int] PIVOTS, negations;
		cdef int _t, t, i, k, j, block, pivot, pivots, M, N, start, stop;
		cdef FFINT q, ratio;

		pivots = 0;
		PIVOTS = Vector[int](self.shape[0]);
		negations = Vector[int](self.shape[0]);

		M = self.shape[0];
		N = self.shape[1];

		# Compute REF.
		for i in range(AUGMENT):
			# Find the row with a pivot in this column (if one exists). If it
			# doesn't, continue; otherwise, invert the row.
			pivot = self.PivotRow(i, pivots)
			if pivot < 0: continue

			# Store the column and swap the rows.
			PIVOTS[pivots] = i
			self.SwapRows(pivots, pivot);

			# Invert the pivot row and increment the number of pivots.
			q = self.inverses[self.data[pivots,i]]
			self.MultiplyRow(pivots, q);
			pivots += 1;

			if not self.parallel:
				for k in range(pivots, M):
					ratio = self.negation[self.data[k,i]]
					self.AddRows(pivots-1, k, i, N, ratio)
			else:
				self.recomputeBlockSchema(i, N);
				for _t in range(pivots, M): negations[_t] = self.negation[self.data[_t, i]]

				for block in prange(self.blockSchema.size(), nogil=True, schedule="dynamic", num_threads=self.cores):
					start = self.blockSchema[block][0];
					stop = self.blockSchema[block][1];
					t = pivots;

					while t < M:
						ratio = negations[t];
						self.ParallelAddRows(pivots-1, t, start, stop, ratio);
						t = t+1

		# Compute RREF.
		cdef int l, m, r, column;
		l = pivots-1;

		# Working backwards from the last nonzero row...
		while l > -1:
			# ... check that we're not getting an out-of-bounds error...
			if l > 0 and PIVOTS[l] < 1:
				l -= 1;
				continue;

			# ... then find the column containing the pivot element...
			column = PIVOTS[l];

			# ... and, iterating over the rows *above* row `l`, eliminate the
			# nonzero entries there...
			if not self.parallel:
				for m in range(l):
					ratio = self.negation[self.data[m,column]]
					self.AddRows(l, m, column, N, ratio)
				
			else:
				self.ScanRow(l);
				self.recomputeBlockSchema(column, N);
				for _t in range(l): negations[_t] = self.negation[self.data[_t, column]]

				for j in prange(self.blockSchema.size(), nogil=True, schedule="dynamic", num_threads=self.cores):
					start = self.blockSchema[j][0];
					stop = self.blockSchema[j][1];
					m = 0;

					while m < l:
						ratio = negations[m]
						self.ParallelAddRows(l, m, start, stop, ratio)
						m = m+1

			# ... then decrement the counter.
			l -= 1

		# Set the highest zero row.
		self.zero = pivots
