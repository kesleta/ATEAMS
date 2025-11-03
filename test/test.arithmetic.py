
from ateams.complexes import Cubical
from ateams.models import SwendsenWang
from ateams.statistics import Chain

field = 3
dimension = 2
N = 100

clex = Cubical().fromCorners([10,10,10,10])
SW = SwendsenWang(clex, dimension=dimension, field=field)
M = Chain(SW, steps=N)

for (spins, occupied) in M:
	pass