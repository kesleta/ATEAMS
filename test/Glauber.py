
from ateams.complexes import Cubical
from ateams.models import Glauber
from ateams.statistics import constant, critical
from ateams import Chain
import json


def construct(L):
	# Construct lattice object.
	field = 3
	L = Cubical().fromCorners([L]*4, field=field)

	# Set up Model and Chain.
	T = critical(L.field)
	SW = Glauber(L, dimension=2, temperature=lambda t: -T(t))
	N = 20
	M = Chain(SW, steps=N)

	return M

def chain(M, DESC=""):
	for (spins, occupied) in M.progress(dynamic_ncols=True, desc=DESC):
		pass


if __name__ == "__main__":
	M = construct(10)
	chain(M)

