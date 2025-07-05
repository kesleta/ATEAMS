
from ateams.complexes import Cubical
from ateams.models import Glauber
from ateams.statistics import constant, critical
from ateams import Chain
import json


def construct(L, field, DIM):
	# Construct complex object.
	L = Cubical().fromCorners([L]*DIM, field=field)

	# Set up Model and Chain.
	T = critical(L.field)
	SW = Glauber(L, dimension=DIM//2, temperature=lambda t: -T(t))
	N = 100
	M = Chain(SW, steps=N)

	return M

def chain(M, DESC=""):
	for (spins, occupied) in M.progress(dynamic_ncols=True, desc=DESC):
		pass
	return M._exitcode


if __name__ == "__main__":
	M = construct(100, 2, 2)
	chain(M)

