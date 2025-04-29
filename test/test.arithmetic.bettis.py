
from ateams.arithmetic import Persistence, FINT
from ateams.structures import Lattice
from galois import GF
import numpy as np
import sys
from tqdm import tqdm

F = 3
L = Lattice().fromCorners([3,3], field=F)
P = Persistence(F, L.flattened)

subcomplex = np.array(range(len(L.flattened)))
# subcomplex = subcomplex[:len(subcomplex)-3]
subcomplex = np.concatenate([subcomplex[:11], subcomplex[12:]])

# print(subcomplex)

g = P.ComputeBettiNumbers(subcomplex)
print(g)
