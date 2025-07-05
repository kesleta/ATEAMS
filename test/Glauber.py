
from ateams.complexes import Cubical
from ateams.models import Glauber
from ateams.statistics import constant, critical
from ateams import Chain
import json


def construct(L, DIM):
	# Construct complex object.
	field = 3
	L = Cubical().fromCorners([L]*DIM, field=field)

	# Set up Model and Chain.
	T = critical(L.field)
	SW = Glauber(L, dimension=DIM//2, temperature=lambda t: -T(t))
	N = 1000
	M = Chain(SW, steps=N)

	return M

def chain(M, DESC=""):
	for (spins, occupied) in M.progress(dynamic_ncols=True, desc=DESC):
		pass


if __name__ == "__main__":
	M = construct(10, 4)
	chain(M)

