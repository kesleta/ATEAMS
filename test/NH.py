
from ateams.complexes import Cubical
from ateams.models import Nienhuis
from ateams.statistics import critical
from ateams import Chain
import json
import sys
from pathlib import Path

def construct(L, field):
	# Construct complex object.
	fname = Path(f"./data/cubical.{L}.2.json")
	if not fname.exists():
		fname.parent.mkdir(exist_ok=True, parents=True)
		L = Cubical().fromCorners([L]*2)
		L.toFile(fname)
	else:
		L = Cubical().fromFile(fname)

	# Set up Model and Chain.
	SW = Nienhuis(L, 1, 1, dimension=2, field=field)
	N = 100
	M = Chain(SW, steps=N)

	return M

def chain(M, DESC=""):
	for (spins, faceEnergy, cellEnergy) in M.progress(dynamic_ncols=True, desc=DESC):
		pass
	return M._exitcode


if __name__ == "__main__":
	M = construct(100, 5)
	chain(M)

