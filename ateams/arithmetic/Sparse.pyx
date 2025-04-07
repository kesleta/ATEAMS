
# cython: language_level=3str, initializedcheck=False, c_api_binop_methods=True, nonecheck=False, profile=True, cdivision=True
# cython: binding=True, linetrace=True
# cython: boundscheck=False, wraparound=False
# define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
# distutils: language=c++

from .common cimport FFINT, FLAT, TABLE
from .common import FINT

import numpy as np
cimport numpy as np

from cython.parallel cimport prange
from libcpp cimport bool
from libcpp.unordered_set cimport unordered_set as Set
from libcpp.vector cimport vector as Vector
from libcpp.unordered_map cimport unordered_map as Map


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
		# the minimum block size, and initialize the block schema.
		self.parallel = parallel;
		self.cores = cores;
		self.minBlockSize = minBlockSize;
		self.maxBlockSize = maxBlockSize;
		self.blockSchema = Vector[Vector[int]]();

		# Initialize the `.rows` and `.columns` properties.
		self.columns = Map[int, Set[int]]();
		self.nonzeroColumns = Map[int, Vector[int]]();
		self.nonzeroColumnCounts = Vector[int](self.shape[0]);
		self._initializeColumns();

		# Keep track of the lowest zero row.
		self.zero = -1
	

	cdef int ThreadsRequired(self, int row) noexcept:
		"""
		Compute the number of threads required.
		"""
		return self.cores;


	cdef Vector[Vector[int]] recomputeBlockSchemaFromRow(self, int row) noexcept:
		cdef int threads, columns, blockWidth, nextIndex, blocks, t;
		cdef Vector[int] BLOCK, INDICES;

		# Clear the current block schema and get the required number of threads.
		self.blockSchema.clear();
		threads = self.ThreadsRequired(row);

		# Get the number of nonzero columns, then split things up that way. We want
		# to use as many threads as possible while re-using as few threads as possible,
		# given the user-provided block parameters.
		columns = self.nonzeroColumnCounts[row];
		blockWidth = (columns+threads-1)/threads;
		blocks = (columns+blockWidth-1)/blockWidth;

		INDICES = Vector[int](blocks+1);
		INDICES[0] = 0;

		# Create the index indicators.
		for t in range(1, blocks+1):
			nextIndex = INDICES[t-1]+blockWidth;
			nextIndex = nextIndex if nextIndex <= columns else columns;
			INDICES[t] = nextIndex;

		# Create the blocks from the index indicators.
		for t in range(1, blocks+1):
			BLOCK = Vector[int](2);
			BLOCK[0] = INDICES[t-1];
			BLOCK[1] = INDICES[t];
			self.blockSchema.push_back(BLOCK);

		return self.blockSchema;
	

	cdef void _initializeColumns(self) noexcept:
		"""
		Initializes the `.rows` and `.columns` data structures.
		"""
		cdef int i, j, M, N, t;
		cdef Vector[int] nonzero;

		M = self.shape[0];
		N = self.shape[1];

		for i in range(M):
			nonzero = Vector[int](N);
			t = 0;

			for j in range(N):
				if self.data[i,j] > 0:
					nonzero[t] = j;
					t += 1;
			
			self.nonzeroColumns[i] = nonzero;
			self.nonzeroColumnCounts[i] = t;
	
	
	cdef void SwapRows(self, int i, int j):
		"""
		Swaps rows `i` and `j`.
		"""
		# Swap entries in the database, re-scanning rows first.
		cdef int _column, column;
		cdef FFINT t;
		cdef Set[int] seen = Set[int]();

		self.ScanRow(i);
		self.ScanRow(j);

		# Swap from row i to row j.
		for _column in range(self.nonzeroColumnCounts[i]):
			column = self.nonzeroColumns[i][_column];
			seen.insert(column);

			t = self.data[i,column];
			self.data[i,column] = self.data[j,column];
			self.data[j,column] = t;

		# Swap from row j to row i.
		for _column in range(self.nonzeroColumnCounts[j]):
			column = self.nonzeroColumns[j][_column];
			if seen.contains(column): continue;
			else: seen.insert(column);

			t = self.data[i,column];
			self.data[i,column] = self.data[j,column];
			self.data[j,column] = t;

		# Swap the column data.
		cdef Vector[int] nonzero;
		cdef int count;

		nonzero = self.nonzeroColumns[i];
		self.nonzeroColumns[i] = self.nonzeroColumns[j];
		self.nonzeroColumns[j] = nonzero;

		count = self.nonzeroColumnCounts[i];
		self.nonzeroColumnCounts[i] = self.nonzeroColumnCounts[j];
		self.nonzeroColumnCounts[j] = count;


	cdef void ParallelAddRows(self, int i, int j, int start, int stop, FFINT ratio) noexcept nogil:
		"""
		Adds rows `i` and `j`.
		"""
		# First, perform all the addition, then scan over both rows to delete fill
		# (extra added zeros). WARNING: this operation is *not commutative*! If
		# we're adding row `i` to row `j`, then we only add *nonzero elements*
		# of row `i` to row `j`.
		cdef int _c, c, t;
		cdef FFINT p, q, mult;

		for _c in range(start, stop):
			c = self.nonzeroColumns[i][_c]

			p = self.data[i,c]
			mult = self.multiplication[p, ratio]
			self.data[j,c] = self.addition[mult,self.data[j,c]]


	cdef void AddRows(self, int i, int j, FFINT ratio) noexcept:
		cdef int _c, column;
		cdef FFINT p, q, mult;

		for _c in range(self.nonzeroColumnCounts[i]):
			column = self.nonzeroColumns[i][_c];
			
			p = self.data[i, column];
			mult = self.multiplication[p,ratio];
			self.data[j,column] = self.addition[mult,self.data[j,column]]
		

	cdef void MultiplyRow(self, int i, FFINT q) noexcept:
		"""
		Multiplies row `i` by `q`.
		"""
		cdef int _c, c;
		cdef FFINT p;

		for _c in range(self.nonzeroColumnCounts[i]):
			c = self.nonzeroColumns[i][_c];
			p = self.data[i,c];
			self.data[i,c] = self.multiplication[q,p];

	
	cdef void ScanRow(self, int i) noexcept:
		cdef int c, t = 0;
		
		for c in range(self.shape[1]):
			if self.data[i,c] > 0:
				self.nonzeroColumns[i][t] = c;
				t += 1

		self.nonzeroColumnCounts[i] = t;

	
	cdef void ScanRowFromColumn(self, int i, int MINCOL) noexcept:
		cdef int c, t = 0;
		
		for c in range(MINCOL, self.shape[1]):
			if self.data[i,c] > 0:
				self.nonzeroColumns[i][t] = c;
				t += 1

		self.nonzeroColumnCounts[i] = t;


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
					p = i if p < 0 and i > pivots-1 else p
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
		cdef int _t, t, i, k, j, threads, block, pivot, pivots, M, N, start, stop;
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
					self.AddRows(pivots-1, k, ratio)
			else:
				self.recomputeBlockSchemaFromRow(pivots-1);
				threads = self.ThreadsRequired(pivots-1);

				for _t in range(pivots, M): negations[_t] = self.negation[self.data[_t, i]]

				for block in prange(self.blockSchema.size(), nogil=True, schedule="dynamic", num_threads=threads):
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
			self.ScanRowFromColumn(l, column);

			# ... and, iterating over the rows *above* row `l`, eliminate the
			# nonzero entries there...
			if not self.parallel:
				for m in range(l):
					ratio = self.negation[self.data[m,column]]
					self.AddRows(l, m, ratio)
				
			else:
				self.recomputeBlockSchemaFromRow(l);
				threads = self.ThreadsRequired(l);

				for _t in range(l): negations[_t] = self.negation[self.data[_t, column]]

				for j in prange(self.blockSchema.size(), nogil=True, schedule="dynamic", num_threads=threads):
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



cdef class ReducedMatrix:
	def __cinit__(
			self,
			TABLE A,
			int characteristic,
			bool parallel,
			int minBlockSize,
			int maxBlockSize,
			int cores
		):
		# Initialize addition and multiplication tables.
		self.characteristic = characteristic;
		self.__arithmetic();
		self.data = A;

		# Initialize the `.shape` property.
		self.shape = Vector[int]();
		self.shape.push_back(A.shape[0]);
		self.shape.push_back(A.shape[1]);

		# If we want synchronous matrix operations, set the number of cores,
		# the minimum block size, and initialize the block schema.
		self.parallel = parallel;
		self.cores = cores;
		self.minBlockSize = minBlockSize;
		self.maxBlockSize = maxBlockSize;
		self.blockSchema = Vector[Vector[int]]();

		# Initialize the `.rows` and `.columns` properties.
		self.columns = Map[int, Set[int]]();
		self.nonzeroColumns = Map[int, Vector[int]]();
		self.nonzeroColumnCounts = Vector[int](self.shape[0]);
		self._initializeColumns();

		# Keep track of the lowest zero row.
		self.zero = -1

	
	cdef void __arithmetic(self) noexcept:
		# Given a field characteristic, construct addition and multiplication
		# tables.
		cdef int p = self.characteristic;
		cdef TABLE addition, multiplication;
		cdef FLAT negation, inverse;
		cdef int i, j;

		addition = np.zeros((p,p), dtype=FINT);
		multiplication = np.zeros((p,p), dtype=FINT);

		# Addition, multiplication tables.
		for i in range(p):
			for j in range(p):
				addition[i,j] = <FFINT>((i+j)%p);
				multiplication[i,j] = <FFINT>((i*j)%p);

		self.addition = addition;
		self.multiplication = multiplication;

		# Negations and inverses.
		negation = np.zeros(p, dtype=FINT);
		inverse = np.zeros(p, dtype=FINT);

		negation[0] = 0;
		inverse[0] = 0;

		for i in range(1, p):
			negation[i] = <FFINT>(p-i);
			inverse[i] = <FFINT>(pow(i, p-2)%p);

		self.negation = negation;
		self.inverse = inverse;
	

	cdef int ThreadsRequired(self, int row) noexcept:
		"""
		Compute the number of threads required.
		"""
		return self.cores;


	cdef Vector[Vector[int]] recomputeBlockSchemaFromRow(self, int row) noexcept:
		cdef int threads, columns, blockWidth, nextIndex, blocks, t;
		cdef Vector[int] BLOCK, INDICES;

		# Clear the current block schema and get the required number of threads.
		self.blockSchema.clear();
		threads = self.ThreadsRequired(row);

		# Get the number of nonzero columns, then split things up that way. We want
		# to use as many threads as possible while re-using as few threads as possible,
		# given the user-provided block parameters.
		columns = self.nonzeroColumnCounts[row];
		blockWidth = (columns+threads-1)/threads;
		blocks = (columns+blockWidth-1)/blockWidth;

		INDICES = Vector[int](blocks+1);
		INDICES[0] = 0;

		# Create the index indicators.
		for t in range(1, blocks+1):
			nextIndex = INDICES[t-1]+blockWidth;
			nextIndex = nextIndex if nextIndex <= columns else columns;
			INDICES[t] = nextIndex;

		# Create the blocks from the index indicators.
		for t in range(1, blocks+1):
			BLOCK = Vector[int](2);
			BLOCK[0] = INDICES[t-1];
			BLOCK[1] = INDICES[t];
			self.blockSchema.push_back(BLOCK);

		return self.blockSchema;
	

	cdef void _initializeColumns(self) noexcept:
		"""
		Initializes the `.rows` and `.columns` data structures.
		"""
		cdef int i, j, M, N, t;
		cdef Vector[int] nonzero;

		M = self.shape[0];
		N = self.shape[1];

		for i in range(M):
			nonzero = Vector[int](N);
			t = 0;

			for j in range(N):
				if self.data[i,j] > 0:
					nonzero[t] = j;
					t += 1;
			
			self.nonzeroColumns[i] = nonzero;
			self.nonzeroColumnCounts[i] = t;
	
	
	cdef void SwapRows(self, int i, int j):
		"""
		Swaps rows `i` and `j`.
		"""
		# Swap entries in the database, re-scanning rows first.
		cdef int _column, column;
		cdef FFINT t;
		cdef Set[int] seen = Set[int]();

		self.ScanRow(i);
		self.ScanRow(j);

		# Swap from row i to row j.
		for _column in range(self.nonzeroColumnCounts[i]):
			column = self.nonzeroColumns[i][_column];
			seen.insert(column);

			t = self.data[i,column];
			self.data[i,column] = self.data[j,column];
			self.data[j,column] = t;

		# Swap from row j to row i.
		for _column in range(self.nonzeroColumnCounts[j]):
			column = self.nonzeroColumns[j][_column];
			if seen.contains(column): continue;
			else: seen.insert(column);

			t = self.data[i,column];
			self.data[i,column] = self.data[j,column];
			self.data[j,column] = t;

		# Swap the column data.
		cdef Vector[int] nonzero;
		cdef int count;

		nonzero = self.nonzeroColumns[i];
		self.nonzeroColumns[i] = self.nonzeroColumns[j];
		self.nonzeroColumns[j] = nonzero;

		count = self.nonzeroColumnCounts[i];
		self.nonzeroColumnCounts[i] = self.nonzeroColumnCounts[j];
		self.nonzeroColumnCounts[j] = count;


	cdef void ParallelAddRows(self, int i, int j, int start, int stop, FFINT ratio) noexcept nogil:
		"""
		Adds rows `i` and `j`.
		"""
		# First, perform all the addition, then scan over both rows to delete fill
		# (extra added zeros). WARNING: this operation is *not commutative*! If
		# we're adding row `i` to row `j`, then we only add *nonzero elements*
		# of row `i` to row `j`.
		cdef int _c, c, t;
		cdef FFINT p, q, mult;

		for _c in range(start, stop):
			c = self.nonzeroColumns[i][_c]

			p = self.data[i,c]
			mult = self.multiplication[p, ratio]
			self.data[j,c] = self.addition[mult,self.data[j,c]]


	cdef void AddRows(self, int i, int j, FFINT ratio) noexcept:
		cdef int _c, column;
		cdef FFINT p, q, mult;

		for _c in range(self.nonzeroColumnCounts[i]):
			column = self.nonzeroColumns[i][_c];
			
			p = self.data[i, column];
			mult = self.multiplication[p,ratio];
			self.data[j,column] = self.addition[mult,self.data[j,column]]
		

	cdef void MultiplyRow(self, int i, FFINT q) noexcept:
		"""
		Multiplies row `i` by `q`.
		"""
		cdef int _c, c;
		cdef FFINT p;

		for _c in range(self.nonzeroColumnCounts[i]):
			c = self.nonzeroColumns[i][_c];
			p = self.data[i,c];
			self.data[i,c] = self.multiplication[q,p];

	
	cdef void ScanRow(self, int i) noexcept:
		cdef int c, t = 0;
		
		for c in range(self.shape[1]):
			if self.data[i,c] > 0:
				self.nonzeroColumns[i][t] = c;
				t += 1

		self.nonzeroColumnCounts[i] = t;

	
	cdef void ScanRowFromColumn(self, int i, int MINCOL) noexcept:
		cdef int c, t = 0;
		
		for c in range(MINCOL, self.shape[1]):
			if self.data[i,c] > 0:
				self.nonzeroColumns[i][t] = c;
				t += 1

		self.nonzeroColumnCounts[i] = t;


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
					p = i if p < 0 and i > pivots-1 else p
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
		cdef int _t, t, i, k, j, threads, block, pivot, pivots, M, N, start, stop;
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
			q = self.inverse[self.data[pivots,i]]
			self.MultiplyRow(pivots, q);
			pivots += 1;

			if not self.parallel:
				for k in range(pivots, M):
					ratio = self.negation[self.data[k,i]]
					self.AddRows(pivots-1, k, ratio)
			else:
				self.recomputeBlockSchemaFromRow(pivots-1);
				threads = self.ThreadsRequired(pivots-1);

				for _t in range(pivots, M): negations[_t] = self.negation[self.data[_t, i]]

				for block in prange(self.blockSchema.size(), nogil=True, schedule="dynamic", num_threads=threads):
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
			self.ScanRowFromColumn(l, column);

			# ... and, iterating over the rows *above* row `l`, eliminate the
			# nonzero entries there...
			if not self.parallel:
				for m in range(l):
					ratio = self.negation[self.data[m,column]]
					self.AddRows(l, m, ratio)
				
			else:
				self.recomputeBlockSchemaFromRow(l);
				threads = self.ThreadsRequired(l);

				for _t in range(l): negations[_t] = self.negation[self.data[_t, column]]

				for j in prange(self.blockSchema.size(), nogil=True, schedule="dynamic", num_threads=threads):
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
