
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
	def __cinit__(self, TABLE A, TABLE addition, FLAT negation, TABLE multiplication, FLAT inverses, bool parallel):
		# Initialize addition and multiplication tables.
		self.addition = addition;
		self.negation = negation;
		self.multiplication = multiplication;
		self.inverses = inverses;
		self.data = A
		self.parallel = parallel

		# Initialize the `.shape` property.
		self.shape = Vector[int]();
		self.shape.push_back(A.shape[0]);
		self.shape.push_back(A.shape[1]);

		# Initialize the `.rows` and `.columns` properties.
		self.columns = Map[int, Set[int]]();
		self.iterableColumns = Map[int, Vector[int]]();
		self._initializeColumns();
	

	cdef void _initializeColumns(self) noexcept:
		"""
		Initializes the `.rows` and `.columns` data structures.
		"""
		cdef int i, j, M, N;
		cdef Set[int] columns;
		cdef Vector[int] iterable;

		M = self.shape[0];
		N = self.shape[1];

		for i in range(M):
			columns = Set[int]();
			iterable = Vector[int]();

			for j in range(N):
				if self.data[i,j] > 0:
					columns.insert(j);
					iterable.push_back(j);
			
			self.columns[i] = columns;
			self.iterableColumns[i] = iterable;

		# self.columns[-1] = Set[int]();
		# self.iterableColumns[-1] = Vector[int]();

	
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

		# # Swap column labels.
		# self.columns[-1] = self.columns[i];
		# self.columns[i] = self.columns[j];
		# self.columns[j] = self.columns[-1];

		# self.iterableColumns[-1] = self.iterableColumns[i];
		# self.iterableColumns[i] = self.iterableColumns[j];
		# self.iterableColumns[j] = self.iterableColumns[-1];

		# Swap the columns (just by swapping keys).
		cdef Set[int] colswap = self.columns[i];
		cdef Vector[int] iterable = self.iterableColumns[i];

		self.columns[i] = self.columns[j];
		self.iterableColumns[i] = self.iterableColumns[j];
		self.columns[j] = colswap;
		self.iterableColumns[j] = iterable;

	
	cdef void AddRows(self, int i, int j, FFINT ratio) noexcept nogil:
		"""
		Adds rows `i` and `j`.
		"""
		# First, perform all the addition, then scan over both rows to delete fill
		# (extra added zeros). WARNING: this operation is *not commutative*! If
		# we're adding row `i` to row `j`, then we only add *nonzero elements*
		# of row `i` to row `j`.
		cdef int k, N, c;
		cdef FFINT p, q, mult;

		N = self.iterableColumns[i].size();

		# if self.parallel:
		# 	for k in prange(N, num_threads=2):
		# 		c = self.iterableColumns[i][k]
		# 		p = self.data[i,c]
		# 		mult = self.multiplication[p, ratio]
		# 		self.data[j,c] = self.addition[mult,self.data[j,c]]

		# 		# Throw out anything that results in zero, and add anything that
		# 		# results in something nonzero.
		# 		if self.data[j,c] < 1: self.columns[j].erase(c)
		# 		else: self.columns[j].insert(c)
		# else:
		for k in range(N):
			c = self.iterableColumns[i][k]
			p = self.data[i,c]
			mult = self.multiplication[p, ratio]
			self.data[j,c] = self.addition[mult,self.data[j,c]]

			# Throw out anything that results in zero, and add anything that
			# results in something nonzero.
			if self.data[j,c] < 1: self.columns[j].erase(c)
			else: self.columns[j].insert(c)

		# Re-do the `j`th column.
		self.iterableColumns[j].clear();
		for c in self.columns[j]: self.iterableColumns[j].push_back(c);


	cdef void MultiplyRow(self, int i, FFINT q) noexcept nogil:
		"""
		Multiplies row `i` by `q`.
		"""
		cdef int k, c, N;
		cdef FFINT p;

		N = self.iterableColumns[i].size();

		if self.parallel:
			for k in prange(N):
				c = self.iterableColumns[i][k];
				p = self.data[i,c];
				self.data[i,c] = self.multiplication[q,p];
		else:
			for k in range(N):
				c = self.iterableColumns[i][k];
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

	
	cdef int HighestZeroRow(self, int AUGMENT=-1) noexcept:
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
			c = min(self.iterableColumns[i])
			if c >= AUGMENT: return i

		return -1

	
	cdef TABLE ToArray(self) noexcept: return self.data
	

	cdef void RREF(self, int AUGMENT=-1) noexcept nogil:
		"""
		Computes the RREF of this matrix.
		"""
		# Catch for optional arguments.
		if AUGMENT < 0: AUGMENT = self.shape[1]

		cdef Vector[int] PIVOTS = Vector[int](self.shape[0]);
		cdef int i, k, pivot, pivots;
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
			if self.parallel:
				for k in prange(pivots, self.shape[0]):
					ratio = self.negation[self.data[k,i]]
					self.AddRows(pivots-1, k, ratio)
			else:
				for k in range(pivots, self.shape[0]):
					ratio = self.negation[self.data[k,i]]
					self.AddRows(pivots-1, k, ratio)
			
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
			if self.parallel:
				for m in prange(l):
					ratio = self.negation[self.data[m,column]]
					self.AddRows(l, m, ratio);
			else:
				for m in range(l):
					ratio = self.negation[self.data[m,column]]
					self.AddRows(l, m, ratio);

			# ... then decrement the counter.
			l -= 1
