
# cython: language_level=3str, initializedcheck=False, c_api_binop_methods=True, nonecheck=False, profile=True, cdivision=True
# cython: binding=True, linetrace=True
# cython: boundscheck=False, wraparound=False
# define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
# distutils: language=c++

from .common cimport TABLE, FLAT, FFINT
from .common import FINT

import numpy as np
cimport numpy as np

from cython.operator cimport dereference
from libcpp.unordered_set cimport unordered_set as Set
from libcpp.set cimport set as OrderedSet
from libcpp.unordered_map cimport unordered_map as Map
from libcpp cimport bool
from libc.math cimport pow


cdef class Persistence:
	def __init__(
			self,
			int homology,
			FFINT characteristic,
			list[list[int]] flattened,
		):
		"""
		Args:
			homology (int): Degree of homology observed for homological percolation events.
			characteristic (int): Characteristic of the finite field over which
				we do computations.
			flattened (list): List-of-lists representing the *original* flattened
				boundary matrix for the complex.
			tranches (np.ndarray): `tranches` property of the `Lattice` object.
			dimensions (np.ndarray): The dimension of each cell.
			
		"""
		self.homology = homology;
		self.characteristic = characteristic;
		self.boundary = self.Vectorize(flattened);

		self.cellCount = self._boundary.size();
		self.vertexCount = self._tranches[0][1];
		self.higherCellCount = self._tranches[self.homology+1][1];
		self.low = self._tranches[self.homology][0];
		self.high = self._tranches[self.homology][1];

		# Construct arithmetic operations.
		self.__arithmetic();

		# Load premarked indices into a `Vector` so we can refresh quickly.
		cdef int i;

		self.premarked = Vector[int](self.vertexCount);
		for i in range(self.vertexCount): self.premarked[i] = i;

		# Data structures for storing column information.
		self.columnEntries = Vector[OrderedSet[int]](self.cellCount);
		self.columnEntriesIterable = Vector[Vector[int]](self.cellCount);
		self.columnEntriesCoefficients = Vector[Map[int,FFINT]](self.cellCount);

		self.__flushDataStructures();


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

	
	cdef void __flushDataStructures(self) noexcept:
		# Data structures for storing column information.
		self.columnEntries = Vector[OrderedSet[int]](self.cellCount);
		self.columnEntriesIterable = Vector[Vector[int]](self.cellCount);
		self.columnEntriesCoefficients = Vector[Map[int,FFINT]](self.cellCount);

		self.marked = Set[int]();
		self.marked.insert(self.premarked.begin(), self.premarked.end());
		self.markedIterable = Vector[int](self.cellCount);


	cpdef Vector[Vector[int]] ReindexBoundary(self, INDEXFLAT filtration) noexcept:
		"""
		Re-indexes the boundary matrix according to the given filtration.
		"""
		cdef int t, i, j, filtered, unfiltered, face, N, M, dimension, start, stop, please;
		cdef Vector[int] faces, indices, temp;

		N = self._boundary.size();
		indices = Vector[int](N);

		# Maps the elements of the filtration to their new indices.
		for t in range(N):
			filtered = filtration[t];
			indices[filtered] = t;

		start = self._tranches[self.homology][0];
		stop = self._tranches[self.homology][1];
		please = self._tranches[self.homology+1][1];

		# Swaps elements when necessary, relabels elements when necessary.
		for i in range(N):
			# We're swapping elements of dimension `self.homology`:
			if i >= start and i < stop:
				filtered = filtration[i];
				self.boundary[i] = self._boundary[filtered];
			# Otherwise, we're relabeling indices.
			elif i >= stop and i < please:
				faces = self._boundary[i];
				M = faces.size();

				for j in range(M):
					unfiltered = faces[j];
					filtered = indices[unfiltered];
					self.boundary[i][j] = filtered;
		
		return self.boundary
		

	cdef Vector[Vector[int]] Vectorize(self, list[list[int]] flattened) noexcept:
		"""
		Convert a list of lists into C++ vectors.

		Args:
			flattened (list): Sparse boundary matrix.

		Returns:
			`std::vector` representing the same data.
		"""
		cdef Vector[Vector[int]] outer;
		cdef Vector[int] inner, dimensions;
		cdef int i, j, lower, M, N, dimension;

		# We only want to get boundary matrices for indices of dimension 1 or
		# greater; i.e. we're excluding vertices.
		M = len(flattened);
		outer = Vector[Vector[int]](M);
		dimensions = Vector[int](M);

		for i in range(M):
			# If we're dealing with a non-vertex, we just copy the original
			# boundary matrix over; otherwise, we're dealing with a vertex, which
			# has no boundary.
			try:
				N = len(flattened[i]);
				inner = Vector[int](N);
				for j in range(N): inner[j] = flattened[i][j];

				dimension = <int>(N/2);
			except:
				inner = Vector[int](0);
				dimension = 0;
			
			dimensions[i] = dimension;
			outer[i] = inner;

		# Set this as the current boundary matrix, but also return it (for
		# completion's sake).
		self._boundary = outer;
		self._dimensions = dimensions;

		# Scan over the array of dimensions, checking how many cells of each
		# dimension exist (and where their indices stop/start).
		cdef Vector[Vector[int]] tranches = Vector[Vector[int]]();
		cdef Vector[int] tranche;
		cdef int t, _tranche = 0;

		N = dimensions.size();

		for t in range(1, N):
			if dimensions[t] > dimensions[t-1] or t == N-1:
				tranche = Vector[int](2);
				tranche[0] = _tranche;
				tranche[1] = t if t < N-1 else N;
				tranches.push_back(tranche);

				_tranche = t;

		self._tranches = tranches;
		
		return outer;


	cdef OrderedSet[int] Eliminate(self, int youngest, OrderedSet[int] faces, Map[int,FFINT] &faceCoefficients) noexcept:
		cdef Vector[int] entriesIterable;
		cdef int i, entry, N;
		cdef FFINT _q, inverse, q, entryCoefficient, faceCoefficient, mul, add;

		# Get the coefficient of the pivot entry, then eliminate the row.
		_q = self.columnEntriesCoefficients[youngest][youngest];
		inverse = self.inverse[_q];

		entriesIterable = self.columnEntriesIterable[youngest];
		N = entriesIterable.size();

		# Eliminate.
		for i in range(N):
			# `entry` is the index of the row (column?) we're trying to eliminate.
			entry = entriesIterable[i];
			entryCoefficient = self.columnEntriesCoefficients[youngest][entry];
			mul = self.multiplication[inverse,entryCoefficient];
			q = self.negation[mul];

			# If `add` is 0, we don't care about tracking it anymore; throw it
			# out of the faces (and remove it from the dict, if it were in there
			# in the first place?).
			if faces.contains(entry):
				faceCoefficient = faceCoefficients[entry];
				add = self.addition[faceCoefficient,q];
			else:
				add = q;

			if add < 1:
				faces.erase(entry);
				faceCoefficients.erase(entry);
			else:
				faces.insert(entry);
				faceCoefficients[entry] = add;
			
		return faces


	cdef OrderedSet[int] RemoveUnmarkedCells(self, int cell, OrderedSet[int] faces, Map[int,FFINT] &faceCoefficients) noexcept:
		"""
		Given the latest cell added to the filtration, remove unmarked
		(nonpivot) entries from its boundary.

		Args:
			int (cell): The most recent cell added to the filtration.

		Returns:
			An `OrderedSet` of `int`s corresponding to (indices of!) pivot
			entries in the boundary of cell `cell`.
		"""
		cdef int i, face, parity;

		for i in range(self.boundary[cell].size()):
			face = self.boundary[cell][i];

			# If the face is unmarked, throw it out.
			if self.marked.contains(face):
				faces.insert(face);

				if (i%2) < 1: faceCoefficients[face] = self.negation[1];
				else: faceCoefficients[face] = self.negation[self.characteristic-1];
		
		return faces;

	
	cdef OrderedSet[int] ReducePivotRow(self, int cell, OrderedSet[int] faces, Map[int,FFINT] &faceCoefficients) noexcept:
		"""
		Reduces the pivot row corresponding to cell `cell`.

		Args:
			cell (int): The most recent cell added to the filtration.

		Returns:
			An ordered set of indices corresponding to nonzero entries in the
			of the pivot column.
		"""
		cdef OrderedSet[int] youngestColumnEntries;
		cdef int youngest;
		cdef FFINT _q, q;

		faces = self.RemoveUnmarkedCells(cell, faces, faceCoefficients);

		while not faces.empty():
			# Get the face of `cell` of maximum degree (i.e. the one added latest
			# to the complex).
			youngest = dereference(faces.rbegin());
			youngestColumnEntries = self.columnEntries[youngest];

			# If the column's empty, we've found a pivot, and we're done.
			if youngestColumnEntries.empty(): break;

			# Otherwise, get the coefficient for the youngest face, and eliminate.
			faces = self.Eliminate(youngest, faces, faceCoefficients);
		
		return faces;
			

	cpdef OrderedSet[int] ComputePercolationEvents(self, INDEXFLAT filtration) noexcept:
		"""
		Computes the times of homological percolation events given a filtration
		and a boundary matrix.

		Args:
			filtration (np.array): Order in which cells are added.
			boundary (list): Flattened boundary matrix.

		Returns:
			A `set` of times at which homological percolation occurs.
		"""
		# Flush the set of marked indices, adding premarked ones.
		self.__flushDataStructures();
		events = OrderedSet[int]();

		# Construct the boundary matrix for this filtration; variables for
		# objects.
		cdef Vector[Vector[int]] boundary = self.ReindexBoundary(filtration);
		cdef Vector[int] facesIterable, degree = Vector[int](self.cellCount);
		cdef OrderedSet[int] faces;
		cdef Map[int,FFINT] faceCoefficients;
		cdef int t, j, cell, time, tagged, youngest;

		# TODO: shouldn't have to iterate over vertices
		tagged = 0;

		for t in range(self.vertexCount, self.higherCellCount):
			# Since our filtrations are discrete (i.e. we add exactly one simplex
			# at each time-step), the "degree" of the cell is the same as the time
			# at which it was added.
			cell = t;

			# Create buckets for indices and coefficients; these are created and
			# stored ONCE, but accessed many times.
			faces = OrderedSet[int]();
			faceCoefficients = Map[int,FFINT]();

			# Eliminate.
			faces = self.ReducePivotRow(cell, faces, faceCoefficients);
			
			# If we end up with no faces, we've found a pivot row.
			if faces.empty():
				# Mark the face, but only add it to the marked iterable if it's
				# of the right dimension.
				self.marked.insert(cell);
				degree[cell] = t;
				
				if self._dimensions[cell] == self.homology:
					self.markedIterable[tagged] = cell;
					tagged = tagged + 1;

			# Otherwise, store the row in the appropriate locations.
			else:
				youngest = dereference(faces.rbegin());
				facesIterable = Vector[int]();
				facesIterable.insert(facesIterable.begin(), faces.begin(), faces.end());

				self.columnEntries[youngest] = faces;
				self.columnEntriesIterable[youngest] = facesIterable;
				self.columnEntriesCoefficients[youngest] = faceCoefficients;
				degree[youngest] = t;

		# # Once we're done eliminating, check over what we find.
		cdef OrderedSet[int] unmarked;

		for j in range(tagged):
			cell = self.markedIterable[j];

			unmarked = self.columnEntries[cell];
			if unmarked.empty(): events.insert(degree[cell]);
		
		return events
