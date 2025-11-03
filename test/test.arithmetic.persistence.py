
import numpy as np

from ateams.common import MINT
from ateams.arithmetic import Twist
from ateams.complexes import Cubical

field = 3
dimension = 1

C = Cubical().fromCorners([3,3])
bd, cbd = C.recomputeBoundaryMatrices(dimension)
T = Twist(field, C.matrices.full, C.breaks, len(C.flattened), dimension)

# T.LinearComputeCobasis()

filtration = np.arange(len(C.flattened), dtype=MINT)
# print(filtration)
# filtration[9] = 10
# filtration[10] = 9

low, high = C.breaks[dimension], C.breaks[dimension+1]
e = T.LinearComputePercolationEvents(filtration);
e = np.array(list(e))
e.sort()
e = e[(e >= low) & (e < high)]

T.RankComputePercolationEvents(filtration);
print(C.breaks)
print(e)
