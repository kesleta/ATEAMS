
from ateams.arithmetic import computeGiantCyclePairs, FINT, boundaryMatrix, Persistence
from ateams.structures import Lattice
from math import comb
from functools import partial
import json
import numpy as np
import sys
from tqdm import tqdm


def constructDefaults(L, sparse, parallel, minBlockSize, maxBlockSize, cores):
	LATTICE = Lattice().fromCorners([L]*4, field=3, dimension=3)

	homology = 2
	t = LATTICE.tranches
	p = LATTICE.field.characteristic

	# Premake the "occupied cells" array; change the dimension of the lattice
	# to correspond to the provided dimension.
	rank = comb(len(LATTICE.corners), homology)
	nullity = len(LATTICE.boundary[homology])

	# Construct a "filtration template" of all the indices. At each step, we
	# determine which of the target indices correspond to satisfied `homology`
	# -dimensional cells, then shuffle the order of those indices.
	targetIndices = np.arange(*t[homology], dtype=FINT)
	zeroedTargetIndices = np.arange(len(LATTICE.boundary[homology]), dtype=FINT)

	# Indices of each cell and dimensions.
	dimensions = np.array(sum([[d]*(t[d,1]-t[d,0]) for d in range(len(t)) if d < homology+2],[]), dtype=FINT)

	times = np.array(range(t[1][0], len(dimensions))).astype(FINT)

	# Find the max over the dimensions and specify a "blank" array the max
	# width.
	zeros = np.zeros((2, t[:,1].max()//8), dtype=FINT)

	# Set multiplicative inverses for the field we're working with.
	fieldInverses = np.array([0] + list(LATTICE.field.Range(1, p)**(-1)), dtype=FINT)

	# Create some pre-fab data structures to provide as fixed arguments to
	# the proposal method.
	low = t[0][1]
	premarked = np.array(list(range(low)), dtype=FINT)

	LATTICE.dimension = homology
	coboundary = boundaryMatrix(LATTICE.boundary, homology, LATTICE.field).T
	identity = np.identity(coboundary.shape[1], dtype=FINT)

	# Create arithmetic lookup tables.
	addition = np.zeros((p,p), dtype=FINT)
	for j in range(p): addition[:,j] = (np.arange(p, dtype=FINT)+j)%p

	subtraction = np.zeros((p,p), dtype=FINT)
	for j in range(p): subtraction[:,j] = (np.arange(p, dtype=FINT)-j)%p

	negation = np.array([-q%p for q in range(0, p)], dtype=FINT)
	inverses = np.array([0] + list(LATTICE.field.Range(1, p)**(-1)), dtype=FINT)

	multiplication = np.zeros((p,p), dtype=FINT)
	for j in range(p): multiplication[:,j] = (np.arange(p, dtype=FINT)*j)%p

	powers = np.full(zeros.shape[1], -1, dtype=FINT)
	powers[1::2] = -powers[1::2]
	powers = powers%p

	pivots = np.zeros(coboundary.shape[1], dtype=FINT)
	result = np.zeros(coboundary.shape[1], dtype=FINT)
	store = np.zeros(coboundary.shape[1], dtype=FINT)
	empty = np.empty(shape=(0,0), dtype=FINT)

	P = Persistence(homology, p, premarked.astype(np.int64), t, dimensions.astype(np.int64))
	
	return partial(
		computeGiantCyclePairs,
		times,
		premarked,
		dimensions,
		t.astype(FINT),
		fieldInverses,
		p,
		homology+1,
		t[homology+1][1],
		np.empty((2,0), dtype=FINT),
		zeros,
		np.empty(zeros.shape[1], dtype=FINT),
		addition.astype(FINT),
		subtraction.astype(FINT),
		multiplication.astype(FINT),
		powers.astype(FINT),
		np.zeros(len(LATTICE.flattened), dtype=FINT)
	), P.ComputePercolationEvents


def test(PERSISTENCE, TESTSET):
	# Test bank.
	with tqdm(total=len(TESTSET)) as bar:
		for i, TEST in enumerate(TESTSET):
			flattened = TEST["flattened"]
			filtration = np.array(TEST["filtration"], dtype=FINT)
			essential = PERSISTENCE(filtration.astype(np.int64), flattened);
			bar.update()