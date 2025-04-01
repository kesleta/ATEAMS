
# cython: language_level=3str, initializedcheck=False, c_api_binop_methods=True, nonecheck=False, profile=True, cdivision=True, boundscheck=False, wraparound=False
# cython: linetrace=True
# cython: binding=True
# define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
# distutils: define_macros=CYTHON_TRACE_NOGIL=1
# distutils: language=c++

from .common cimport FFINT, FLAT, TABLE

from cython.parallel cimport prange
from libcpp cimport bool
from libcpp.unordered_set cimport unordered_set as Set
from libcpp.vector cimport vector as Vector
from libcpp.unordered_map cimport unordered_map as Map


cdef Set[int] _intersect = Set[int]();
cdef Set[int] _union = Set[int]();


cdef Set[int] Union(Set[int] P, Set[int] Q) noexcept nogil:
	cdef Set[int] S = Set[int](P);
	cdef int s;

	for s in Q: S.insert(s);

	return S;


cdef Set[int] Intersection(Set[int] P, Set[int] Q) noexcept nogil:
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
			int cores,
			str schedule
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
		self.schedule = schedule;
		self.blockSchema = Vector[Vector[int]]();

		cdef int _q, b, block, blocks, MINCOL, MAXCOL;
		cdef Vector[int] BLOCK;

		_q = self.shape[1]//self.cores;
		block = min(max(self.minBlockSize, _q), min(self.maxBlockSize, _q));
		blocks = self.shape[1]//block;
		self.blockSize = block;

		for b in range(blocks):
			MINCOL = b*block;
			MAXCOL = (b+1)*block if (b+1)*block <= self.shape[1] else self.shape[1];
			
			if MINCOL >= self.shape[1]: break;

			BLOCK = Vector[int]();
			BLOCK.push_back(MINCOL);
			BLOCK.push_back(MAXCOL);
			self.blockSchema.push_back(BLOCK);

		# Initialize the `.rows` and `.columns` properties.
		self.columns = Map[int, Set[int]]();
		self._initializeColumns();
	

	cdef void _initializeColumns(self) :
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

	
	cdef void SwapRows(self, int i, int j) noexcept nogil:
		"""
		Swaps rows `i` and `j`.
		"""
		# Swap entries in the database.
		cdef int c;
		cdef FFINT t;
		cdef Set[int] COLUMNS = Union(self.columns[i], self.columns[j]);

		for c in COLUMNS:
			t = self.data[i,c];
			self.data[i,c] = self.data[j,c];
			self.data[j,c] = t;

		# Swap the columns (just by swapping keys).
		cdef Set[int] colswap = self.columns[i];
		self.columns[i] = self.columns[j];
		self.columns[j] = colswap;

	
	cdef void AddRows(self, int i, int j, int MINCOL, int MAXCOL, FFINT ratio) noexcept nogil:
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
			if c < MINCOL or c >= MAXCOL: continue

			p = self.data[i,c]
			mult = self.multiplication[p, ratio]
			self.data[j,c] = self.addition[mult,self.data[j,c]]

			# If we aren't doing things in parallel, we can write to the column
			# dictionaries. Otherwise, we just re-scan the array and modify the
			# sets where necessary.
			if not self.parallel:
				# Throw out anything that results in zero, and add anything that
				# results in something nonzero.
				if self.data[j,c] < 1: self.columns[j].erase(c)
				else: self.columns[j].insert(c)


	cdef void RescanColumns(self, int MINCOL) noexcept nogil:
		"""
		If we aren't doing things in parallel, we need to re-scan the array.
		"""
		cdef int i, j;

		for i in range(self.shape[0]):
			for j in range(MINCOL, self.shape[1]):
				if self.data[i,j] < 1: self.columns[i].erase(j);
				else: self.columns[i].insert(j);


	cdef void MultiplyRow(self, int i, FFINT q) noexcept nogil:
		"""
		Multiplies row `i` by `q`.
		"""
		cdef int c;
		cdef FFINT p;

		for c in self.columns[i]:
			p = self.data[i,c];
			self.data[i,c] = self.multiplication[q,p];
		

	cdef int PivotRow(self, int c, int pivots) noexcept nogil:
		"""
		Report the first nonzero entry in column `c`; if no such entry exists,
		return -1.
		"""
		cdef int i;

		for i in range(self.shape[0]):
			if self.columns[i].contains(c) and i > pivots-1: return i

		return -1

	
	cdef int HighestZeroRow(self, int AUGMENT=-1) :
		"""
		Report the first lowest-indexed zero row; if no such row exists,
		report -1. If nonnegative AUGMENT is passed, we check for the first row
		with AUGMENT leading zeros.
		"""
		# Catch for optional arguments.
		if AUGMENT < 0: AUGMENT = self.shape[1]
		cdef int i, c;
		cdef bool single;

		for i in range(self.shape[0]):
			c = min(self.columns[i])
			if c >= AUGMENT: return i

		return -1


	cdef void BlockAddRows(
		self,
		int i,
		int j,
		FFINT ratio,
		int start,
		int stop,
		str schedule
	) noexcept nogil:
		cdef int k, N, first, MINCOL, MAXCOL;

		# Determine the lowest index required by `start` based on the block size.
		first = start//self.blockSize;
		N = self.blockSchema.size();

		# Compute block-wise.
		for k in prange(first, N, schedule="static", nogil=True, num_threads=self.cores, chunksize=1):
			MINCOL = self.blockSchema[k][0]
			MAXCOL = self.blockSchema[k][1]

			self.AddRows(i, j, MINCOL, MAXCOL, ratio);

		self.RescanColumns(start);

	
	cdef TABLE ToArray(self) : return self.data
	

	cdef void RREF(self, int AUGMENT=-1) noexcept nogil:
		"""
		Computes the RREF of this matrix.
		"""
		# Catch for optional arguments.
		if AUGMENT < 0: AUGMENT = self.shape[1]

		cdef Vector[int] PIVOTS = Vector[int](self.shape[0]);
		cdef int i, k, j, pivot, pivots;
		cdef FFINT q, ratio;

		pivots = 0;

		# Compute REF.
		for i in range(AUGMENT):
			# Find the row with a pivot in this column (if one exists). If it
			# doesn't, continue; otherwise, invert the row.
			pivot = self.PivotRow(i, pivots)
			if pivot < 0: continue

			# Store the column and swap row entries.
			PIVOTS[pivots] = i
			self.SwapRows(pivots, pivot);

			# Invert the pivot row and increment the number of pivots.
			q = self.inverses[self.data[pivots,i]]
			self.MultiplyRow(pivots, q);
			pivots += 1;

			# Eliminate rows below.
			for k in range(pivots, self.shape[0]):
				ratio = self.negation[self.data[k,i]]

				# self.AddRows(pivots-1, k, 0, self.shape[1], ratio)
				if not self.parallel: self.AddRows(pivots-1, k, i, self.shape[1], ratio)
				else: self.BlockAddRows(pivots-1, k, ratio, i, self.shape[1], self.schedule)

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
			for m in range(l):
				ratio = self.negation[self.data[m,column]]

				# self.AddRows(l, m, 0, self.shape[1], ratio)
				if not self.parallel: self.AddRows(l, m, column, self.shape[1], ratio)
				else: self.BlockAddRows(l, m, ratio, column, self.shape[1], self.schedule)

			# ... then decrement the counter.
			l -= 1
