

import numpy as np
cimport numpy as np
from ..common cimport FFINT, FLAT, TABLE, FLATCONTIG, TABLECONTIG
from ..common import FINT
from .MatrixReduction cimport MatrixReduction


cpdef TABLECONTIG Kernel(MatrixReduction M, TABLE A) noexcept
cpdef np.ndarray[FFINT, ndim=1, negative_indices=False, mode="c"] KernelSample(MatrixReduction M, TABLE A) noexcept
