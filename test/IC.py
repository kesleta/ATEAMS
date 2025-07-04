
from ateams.complexes import Cubical
from ateams.models import InvadedCluster
from ateams import Chain
import json
import sys


def construct(L, parallel, cores, LinBox):
	# Construct lattice object.
	field = 3
	L = Cubical().fromCorners([L]*4, field=field)

	# Set up Model and Chain.
	homology = 2
	SW = InvadedCluster(L, dimension=homology, parallel=parallel, cores=cores, LinBox=LinBox)
	N = 10
	M = Chain(SW, steps=N)

	return M

def chain(M, DESC=""):
	try:
		for (spins, essentials, satisfied) in M.progress(dynamic_ncols=True, desc=DESC):
			pass
	except Exception as e:
		print(e)
		exit(1)

if __name__ == "__main__":
	M = construct(3, False, 2, True)
	chain(M)

