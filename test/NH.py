
from ateams.complexes import Cubical
from ateams.models import Nienhuis
from ateams.statistics import critical
from ateams import Chain
import json
import sys

def construct(L, field):
	# Construct complex object.
	L = Cubical().fromCorners([L]*2, field=field)

	# Set up Model and Chain.
	SW = Nienhuis(L, 1, 1, dimension=2)
	N = 50
	M = Chain(SW, steps=N)

	return M

def chain(M, DESC=""):
	for (spins, faceEnergy, cellEnergy) in M.progress(dynamic_ncols=True, desc=DESC):
		pass


if __name__ == "__main__":
	M = construct(100, False, 2, True)
	chain(M)

