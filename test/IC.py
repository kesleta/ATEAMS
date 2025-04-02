
from ateams.structures import Lattice
from ateams.models import CInvadedCluster, InvadedCluster
from ateams import Chain


def construct(L, sparse, parallel, minBlockSize, maxBlockSize, cores):
	# Construct lattice object.
	field = 3
	L = Lattice().fromCorners([L]*4, dimension=3, field=field)

	# Set up Model and Chain.
	homology = 2
	SW = CInvadedCluster(L, homology=homology, sparse=sparse, parallel=parallel, minBlockSize=minBlockSize, maxBlockSize=maxBlockSize, cores=cores)
	N = 3
	M = Chain(SW, steps=N)

	return M

def chain(M):
	for (spins, essentials, satisfied) in M.progress():
		pass

if __name__ == "__main__":
	chain()

