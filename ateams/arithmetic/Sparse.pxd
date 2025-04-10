
from .common cimport FFINT, FLAT, TABLE, FLATCONTIG, TABLECONTIG

from libcpp cimport bool
from libcpp.vector cimport vector as Vector
from libcpp.unordered_set cimport unordered_set as Set
from libcpp.unordered_map cimport unordered_map as Map


cdef class MatrixReduction:
	cdef TABLECONTIG addition
	cdef FLATCONTIG negation
	cdef TABLECONTIG multiplication
	cdef FLATCONTIG inverse
	cdef TABLECONTIG data
	cdef bool parallel
	cdef int characteristic
	cdef int zero
	cdef int cores
	cdef int minBlockSize
	cdef int maxBlockSize
	cdef int blockSize

	cdef void __arithmetic(self) noexcept
	
	cdef Map[int, Vector[int]] nonzeroColumns
	cdef Vector[int] nonzeroColumnCounts
	cdef Map[int, Set[int]] columns
	cdef Vector[int] shape
	cdef Vector[Vector[int]] blockSchema

	cdef void _initializeColumns(self, TABLECONTIG A) noexcept
	cpdef TABLECONTIG ToArray(self) noexcept

	cdef void SwapRows(self, int i, int j)
	cdef void AddRows(self, int i, int j, FFINT ratio) noexcept
	cdef void ParallelAddRows(self, int i, int j, int start, int stop, FFINT ratio) noexcept nogil
	cdef void MultiplyRow(self, int i, FFINT q) noexcept
	cdef void ScanRow(self, int i) noexcept
	cdef void ScanRowFromColumn(self, int i, int MINCOL) noexcept
	cdef Vector[Vector[int]] recomputeBlockSchemaFromRow(self, int row) noexcept
	cdef int ThreadsRequired(self, int row) noexcept
	cdef int PivotRow(self, int c, int pivots) noexcept
	
	cpdef int HighestZeroRow(self, int AUGMENT=*) noexcept
	cpdef TABLECONTIG RREF(self, TABLECONTIG A, int AUGMENT=*) noexcept


# cdef class Matrix:
# 	cdef TABLE addition
# 	cdef FLAT negation
# 	cdef TABLE multiplication
# 	cdef FLAT inverses
# 	cdef TABLE data
# 	cdef bool parallel
# 	cdef int zero
# 	cdef int cores
# 	cdef int minBlockSize
# 	cdef int maxBlockSize
# 	cdef int blockSize
	
# 	cdef Map[int, Vector[int]] nonzeroColumns
# 	cdef Vector[int] nonzeroColumnCounts
# 	cdef Map[int, Set[int]] columns
# 	cdef Vector[int] shape
# 	cdef Vector[Vector[int]] blockSchema

# 	cdef void _initializeColumns(self) noexcept
# 	cdef TABLE ToArray(self) noexcept

# 	cdef void SwapRows(self, int i, int j)
# 	cdef void AddRows(self, int i, int j, FFINT ratio) noexcept
# 	cdef void ParallelAddRows(self, int i, int j, int start, int stop, FFINT ratio) noexcept nogil
# 	cdef void MultiplyRow(self, int i, FFINT q) noexcept
# 	cdef void ScanRow(self, int i) noexcept
# 	cdef void ScanRowFromColumn(self, int i, int MINCOL) noexcept
# 	cdef Vector[Vector[int]] recomputeBlockSchemaFromRow(self, int row) noexcept
# 	cdef int ThreadsRequired(self, int row) noexcept
# 	cdef int PivotRow(self, int c, int pivots) noexcept
# 	cdef int HighestZeroRow(self, int AUGMENT=*) noexcept

# 	cdef void RREF(self, int AUGMENT=*) noexcept

