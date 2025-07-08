
from ateams.complexes import Cubical
from ateams.models import InvadedCluster
from ateams import Chain
import json
import sys
from pathlib import Path


def construct(L, dim, field):
	# Construct complex object.
	fname = Path(f"./data/cubical.{L}.{dim}.json")
	if not fname.exists():
		fname.parent.mkdir(exist_ok=True, parents=True)
		L = Cubical().fromCorners([L]*dim)
		L.toFile(fname)
	else:
		L = Cubical().fromFile(fname)

	# Set up Model and Chain.
	SW = InvadedCluster(L, dimension=dim//2, field=field)
	N = 20
	M = Chain(SW, steps=N)

	return M

def chain(M, DESC=""):
	for result in M.progress(dynamic_ncols=True, desc=DESC):
		pass
	return M._exitcode

if __name__ == "__main__":
	M = construct(6, 4, 3)
	chain(M)

