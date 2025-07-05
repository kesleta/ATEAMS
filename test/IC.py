
from ateams.complexes import Cubical
from ateams.models import InvadedCluster
from ateams import Chain
import json
import sys


def construct(L, dim, field):
	# Construct complex object.
	L = Cubical().fromCorners([L]*dim, field=field)

	# Set up Model and Chain.
	SW = InvadedCluster(L, dimension=dim//2)
	N = 50
	M = Chain(SW, steps=N)

	return M

def chain(M, DESC=""):
	for (spins, essentials, satisfied) in M.progress(dynamic_ncols=True, desc=DESC):
		pass

if __name__ == "__main__":
	M = construct(5, False, 2, True)
	chain(M)

