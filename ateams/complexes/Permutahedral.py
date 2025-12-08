import json
import math as math
import numpy as np
from itertools import combinations, product, accumulate
from collections import defaultdict

from .construction import fullBoundaryMatrix
from ..common import MINT, Matrices

def dbg(x):
	print(x)
	return x

class Permutahedral:
	def __init__(self): pass

	@classmethod
	def permutahedron(cls, D):
		boundries = {d: {} for d in range(D + 1)}
		boundries[D][(0,) * (D + 1)] = []

		for d in reversed(range(1, D + 1)):
			for cell in boundries[d]:
				for part in range(D + 1 - d):
					indices = [j for j, x in enumerate(cell) if x == part]
					for k in range(1,len(indices)):
						for comb in combinations(indices, k):
							cand = tuple(x+1 if x > part or j in comb else x for j, x in enumerate(cell))
							boundries[d][cell].append(cand)
							boundries[d-1][cand] = []
		return boundries

	def construct(self, dimensions, alternating = False, periodic = True):
		self.D = len(dimensions)
		self.dimensions = dimensions
		self.corners = self.dimensions
		perm = self.permutahedron(self.D)
		centroids = {}

		factor = 2 * math.factorial(self.D + 1)
		inds = [(factor // 2) * i for i in range(-self.D, (self.D+1), 2)]

		for d in range(self.D + 1):
			for cell in perm[d]:
				blocks = [0,] * (self.D - d + 1)
				for p in cell:
					blocks[p] += 1
				starts = [0] + [int(x) for x in np.cumsum(blocks)]
				avgs = [
					sum(inds[starts[i]:starts[i + 1]]) // (starts[i + 1] - starts[i])
					for i in range(self.D - d)
				] + [
					sum(inds[starts[self.D - d]:]) // (len(inds) - starts[self.D - d])
				]
				centroids[cell] = [avgs[p] for p in cell]

		# ###

		# coords = np.array([
		# 	np.array([factor * ((D+1)*lat_coords[i] - sum(lat_coords)) for i in range(D)] + [-factor * sum(lat_coords), i])
		# 	for i, lat_coords in enumerate(product(*(range(d) for d in dimensions)))
		# ])

		# print(coords)

		# ###

		cells_at = defaultdict(lambda: defaultdict(list))

		for lat_coords in product(*(range(d) for d in dimensions)):
			if alternating:
				raise Exception("Haven't implemented alternating yet")
			else:
				coords = tuple(factor * ((self.D+1)*lat_coords[i] - sum(lat_coords)) for i in range(self.D)) + (-factor * sum(lat_coords),)
			for d in range(self.D + 1):
				for cell in perm[d]:
					pos = tuple(x + y for (x, y) in zip(coords, centroids[cell]))
					cells_at[d][pos].append((coords, cell))

		if periodic:
			raise Exception("Haven't implemented perodic yet")

		points = {d: sorted(list(cells_at[d])) for d in range(self.D + 1)}
		point_indices = {d: {} for d  in range(self.D + 1)}

		for d in range(self.D + 1):
			for i, point in enumerate(points[d]):
				point_indices[d][point] = i

		self.Boundary = {d: [[] for _ in points[d]] for d in range(self.D + 1)}

		for d in reversed(range(0, self.D + 1)):
			for point in cells_at[d]:
				(coords, cell) = cells_at[d][point][0]
				self.Boundary[d][point_indices[d][point]] = sorted([
					point_indices[d-1][tuple(x + y for (x, y) in zip(coords, centroids[b]))]
					for b in perm[d][cell]
				])
		
		self.matrices = Matrices()
		self.matrices.boundary, self.matrices.coboundary = self.recomputeBoundaryMatrices(self.D)

		self.breaks = np.array([0] + list(accumulate(len(self.Boundary[d]) for d in range(self.D))), dtype = MINT)
		self.flattened = [(int(c+self.breaks[d]), int(r+self.breaks[d-1])) for d in range(self.D+1) for (c, v) in enumerate(self.Boundary[d]) for r in v]

		coefficients = {
			d: [(-1)**j for j in range(len(self.Boundary[d][0]))]
			for d in range(self.D + 1)
		}
		self.matrices.full = np.array(fullBoundaryMatrix(self.flattened, coefficients), dtype=MINT)
		
		return self


	def recomputeBoundaryMatrices(self, dimension):
		mat = [[c, r, (-1) ** p] for (c, v) in enumerate(self.Boundary[dimension]) for p, r in enumerate(v)]
		comat = [[r, c, x] for [c, r, x] in mat]
		return np.array(sum(mat, []), dtype = MINT), np.array(sum(comat, []), dtype = MINT)
