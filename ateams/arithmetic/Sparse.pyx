
# cython: language_level=3str, initializedcheck=False, c_api_binop_methods=True, nonecheck=False, profile=True, cdivision=True, wraparound=False, boundscheck=False
# define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
# distutils: define_macros=CYTHON_TRACE_NOGIL=1
# distutils: language=c++

from .common cimport FFINT, FLAT, TABLE

from libcpp.set cimport set as Set
from libcpp.vector cimport vector as Vector
from libcpp.map cimport map as Map


cdef Set[int] Intersection(Set[int] P, Set[int] Q):
	cdef Set[int] C, T, S;
	cdef int s;

	C = Set[int]();

	# Only iterate over the smaller set.
	if P.size() < Q.size():
		T = P;
		S = Q;
	else:
		T = Q;
		S = P;

	for s in T:
		if S.contains(s): C.insert(s)

	return C


cdef class Matrix:
	def __cinit__(self, TABLE A, TABLE addition, FLAT negation, TABLE multiplication, FLAT inverses):
		# Initialize addition and multiplication tables.
		self.addition = addition;
		self.negation = negation;
		self.multiplication = multiplication;
		self.inverses = inverses;
		self._original = A;

		# Initialize the `.shape` property.
		self.shape = Vector[int]();
		self.shape.push_back(A.shape[0]);
		self.shape.push_back(A.shape[1]);

		# Initialize the `.rows` and `.columns` properties.
		self.rows = Map[int, Map[int, FFINT]]();
		self.columns = Map[int, Set[int]]();
		self._initializeRows(A);
	

	cdef void _initializeRows(self, TABLE A) noexcept:
		"""
		Initializes the `.rows` and `.columns` data structures.
		"""
		cdef int i, j, M, N;
		cdef Set[int] columns;
		cdef Map[int, FFINT] row;

		M = A.shape[0];
		N = A.shape[1];

		for i in range(M):
			row = Map[int, FFINT]();
			columns = Set[int]();

			for j in range(N):
				if A[i,j] > 0:
					columns.insert(j)
					row[j] = A[i,j]

			self.rows[i] = row
			self.columns[i] = columns

	
	cdef void SwapRows(self, int i, int j) noexcept nogil:
		"""
		Swaps rows `i` and `j`.
		"""
		cdef Map[int, FFINT] rowswap = self.rows[i];
		cdef Set[int] colswap = self.columns[i];

		self.rows[i] = self.rows[j];
		self.rows[j] = rowswap;

		self.columns[i] = self.columns[j];
		self.columns[j] = colswap;

	
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

		# N = self.columns[i].size();

		for c in self.columns[i]:
			p = self.rows[i][c]
			mult = self.multiplication[p, ratio]
			self.rows[j][c] = self.addition[mult,self.rows[j][c]]

			# Throw out anything that results in zero, and add anything that
			# results in something nonzero.
			if self.rows[j][c] < 1:
				self.rows[j].erase(c)
				self.columns[j].erase(c)
			else:
				self.columns[j].insert(c)


	cdef void MultiplyRow(self, int i, FFINT q) noexcept nogil:
		"""
		Multiplies row `i` by `q`.
		"""
		cdef int c;
		cdef FFINT p;

		for c in self.columns[i]:
			p = self.rows[i][c]
			self.rows[i][c] = self.multiplication[q,p]
		

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
		cdef Set[int] COLUMNS = Set[int](range(AUGMENT));

		for i in range(self.shape[0]):
			if Intersection(self.columns[i], COLUMNS).empty(): return i

		return -1

	
	cdef TABLE ToArray(self) noexcept nogil:
		# Overwrites the contents of the original matrix with whatever is stored
		# in the current data table.
		cdef int r, a;
		cdef FFINT p;

		for r in range(self.shape[0]):
			for c in range(self.shape[1]):
				if self.columns[r].contains(c): p = self.rows[r][c]
				else: p = 0

				self._original[r,c] = p

		return self._original
	

	cdef void RREF(self, int AUGMENT=-1) noexcept:
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
			q = self.inverses[self.rows[pivots][i]]
			self.MultiplyRow(pivots, q);
			pivots += 1;

			# Eliminate rows below.
			for k in range(pivots, self.shape[0]):
				ratio = self.negation[self.rows[k][i]]
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
			for m in range(l):
				ratio = self.negation[self.rows[m][column]]
				self.AddRows(l, m, ratio);

			# ... then decrement the counter.
			l -= 1
