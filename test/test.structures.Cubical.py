
from ateams.complexes import Cubical

L = Cubical().fromCorners([3,3,3,3], dimension=2, field=3)
print(L.matrices.boundary)

L.recomputeBoundaryMatrices(1)
print(L.matrices.boundary)