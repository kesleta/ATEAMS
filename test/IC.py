
from ateams.structures import Lattice
from ateams.models import CInvadedCluster, InvadedCluster, CPInvadedCluster
from ateams import Chain
import json


def construct(L, sparse, parallel, minBlockSize, maxBlockSize, cores):
	# Construct lattice object.
	field = 3
	L = Lattice().fromCorners([L]*4, dimension=3, field=field)

	# Set up Model and Chain.
	homology = 2
	SW = CPInvadedCluster(L, homology=homology, sparse=sparse, parallel=parallel, minBlockSize=minBlockSize, maxBlockSize=maxBlockSize, cores=cores)
	N = 20
	M = Chain(SW, steps=N)

	return M

def chain(M):
	for (spins, essentials, satisfied) in M.progress():
		pass

if __name__ == "__main__":
	M = construct(3, False, False, 32, 64, 4)
	chain(M)

