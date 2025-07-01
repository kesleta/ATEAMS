
from ateams.complex import Lattice
from ateams.model import SwendsenWang
from ateams.stats import constant, critical
from ateams import Chain
import json


def construct(L, sparse, parallel, minBlockSize, maxBlockSize, cores, LinBox):
	# Construct lattice object.
	field = 3
	L = Lattice().fromCorners([L]*4, dimension=2, field=field)

	# Set up Model and Chain.
	T = critical(L.field.characteristic)
	SW = SwendsenWang(L, temperatureFunction=lambda t: -T(t), sparse=sparse, parallel=parallel, minBlockSize=minBlockSize, maxBlockSize=maxBlockSize, cores=cores, LinBox=LinBox)
	N = 20
	M = Chain(SW, steps=N)

	return M

def chain(M, DESC=""):
	for (spins, occupied) in M.progress(dynamic_ncols=True, desc=DESC):
		pass

if __name__ == "__main__":
	M = construct(5, False, True, 32, 64, 2, True)
	chain(M)

