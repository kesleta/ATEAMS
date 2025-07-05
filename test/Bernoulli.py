
from ateams.complexes import Cubical
from ateams.models import Bernoulli
from ateams import Chain


def construct(L, parallel, cores, LinBox):
	# Construct complex object.
	field = 2
	L = Cubical().fromCorners([L]*4, field=field)

	# Set up Model and Chain.
	homology = 2
	SW = Bernoulli(L, dimension=homology, parallel=parallel, cores=cores, LinBox=LinBox)
	N = 50
	M = Chain(SW, steps=N)

	return M

def chain(M, DESC=""):
	for essentials in M.progress(dynamic_ncols=True, desc=DESC):
		pass

if __name__ == "__main__":
	M = construct(4, False, 2, True)
	chain(M)

