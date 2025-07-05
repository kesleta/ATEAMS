
from ateams.complexes import Cubical
from ateams.models import Bernoulli
from ateams import Chain


def construct(L, DIM):
	# Construct complex object.
	L = Cubical().fromCorners([L]*DIM)

	# Set up Model and Chain.
	SW = Bernoulli(L, dimension=DIM//2)
	N = 50
	M = Chain(SW, steps=N)

	return M

def chain(M, DESC=""):
	for essentials in M.progress(dynamic_ncols=True, desc=DESC):
		pass

if __name__ == "__main__":
	M = construct(4, 4)
	chain(M)

