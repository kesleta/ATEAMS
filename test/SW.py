
from ateams.structures import Lattice
from ateams.models import SwendsenWang
from ateams.stats import constant, critical
from ateams import Chain
import json


def construct(L, sparse, parallel, minBlockSize, maxBlockSize, cores):
	# Construct lattice object.
	field = 3
	L = Lattice().fromCorners([64]*2, dimension=1, field=field)

	# Set up Model and Chain.
	T = critical(L.field.characteristic)
	SW = SwendsenWang(L, temperatureFunction=lambda t: -T(t))
	N = 20
	M = Chain(SW, steps=N)

	return M

def chain(M, DESC=""):
	for (spins, occupied) in M.progress(dynamic_ncols=True, desc=DESC):
		pass

if __name__ == "__main__":
	M = construct(3, False, True, 32, 64, 2)
	chain(M)

