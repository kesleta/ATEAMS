
import numpy as np
from ateams.complexes import Cubical
from ateams.models import InvadedCluster
from ateams import Chain

L = Cubical().fromCorners([5,5,5,5], field=3)

SW = InvadedCluster(L, dimension=2, LinBox=True)
N = 10
M = Chain(SW, steps=N)

for (spins, occupied, satisfied) in M.progress():
	pass
