
from ateams.complexes import Cubical
from ateams.models import SwendsenWang
from ateams.statistics import constant, critical
from ateams import Chain
import json


def construct(L, dim, field):
	# Construct complex object.
	L = Cubical().fromCorners([L]*dim, field=field)

	# Set up Model and Chain.
	T = critical(L.field)
	SW = SwendsenWang(L, dimension=dim//2, temperature=lambda t: -T(t))
	N = 50
	M = Chain(SW, steps=N)

	return M

def chain(M, DESC=""):
	for (spins, occupied) in M.progress(dynamic_ncols=True, desc=DESC):
		pass


if __name__ == "__main__":
	M = construct(3, 4, 7)
	chain(M)

