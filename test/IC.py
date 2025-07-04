
from ateams.complexes import Cubical
from ateams.models import InvadedCluster
from ateams import Chain
import json


def construct(L, sparse, parallel, minBlockSize, maxBlockSize, cores):
	# Construct lattice object.
	field = 3
	L = Cubical().fromCorners([L]*4, field=field)

	# Set up Model and Chain.
	homology = 2
	SW = InvadedCluster(L, dimension=homology)
	N = 20
	M = Chain(SW, steps=N)

	return M

def chain(M, DESC=""):
	for (spins, essentials, satisfied) in M.progress(dynamic_ncols=True, desc=DESC):
		pass

if __name__ == "__main__":
	M = construct(3, False, False, 32, 64, 2)
	chain(M)

