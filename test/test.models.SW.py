
from ateams.complexes import Cubical
from ateams.models import SwendsenWang
from ateams.statistics import constant
from ateams import Chain

L = Cubical().fromCorners([8]*4, field=3)

SW = SwendsenWang(L, dimension=2, temperature=constant(-0.585))
N = 1000
M = Chain(SW, steps=N)

for (spins, occupied) in M.progress():
	pass
