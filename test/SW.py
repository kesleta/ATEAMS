
from ateams.complexes import Cubical
from ateams.models import SwendsenWang
from ateams.statistics import constant, critical
from ateams import Chain
import json


def construct(L, parallel, cores, LinBox):
	# Construct complex object.
	field = 3
	L = Cubical().fromCorners([L]*4, field=field)

	# Set up Model and Chain.
	T = critical(L.field)
	SW = SwendsenWang(L, dimension=2, temperature=lambda t: -T(t), LinBox=LinBox, parallel=parallel, cores=cores)
	N = 100
	M = Chain(SW, steps=N)

	return M

def chain(M, DESC=""):
	for (spins, occupied) in M.progress(dynamic_ncols=True, desc=DESC):
		pass


if __name__ == "__main__":
	M = construct(7, False, 2, True)
	chain(M)

