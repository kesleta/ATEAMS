
import numpy as np

from ateams.common import MINT
from ateams.arithmetic import ComputePercolationEvents
from ateams.complexes import Cubical

C = Cubical().fromCorners([3,3], field=3)
_, __, full = C.recomputeBoundaryMatrices(2)

filtration = np.arange(len(C.flattened), dtype=MINT)
filtration[9] = 10
filtration[10] = 9
print(C.flattened)

s = ComputePercolationEvents(full, filtration, 1, C.field, C.breaks)
print(s)
