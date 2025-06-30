
from ateams.arithmetic import Fast, FINT
import numpy as np
from galois import GF

p = 3
F = GF(p)

shapes = [(50,60), (60,50), (50,50)]

A = F.Random(shapes[1])
# A = [
# 	[0,0],
# 	[0,1]
# ]
A = np.array(A, dtype=FINT)

g = Fast(A, p)
print(g)
