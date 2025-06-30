
from .common import FINT, SFINT, MINT
from .Persistence import Persistence
from .MatrixReduction import MatrixReduction
from .linalg import Kernel, KernelSample
from .cubicalComplex import cubicalComplex, boundaryMatrix, flatten
from .fastiteration import energy
from .reindexing import reindexSparseBoundaryMatrix
from .linearAlgebra import autocorrelation, evaluateCochain, isNullHomologous
from .Fast import Fast


__pdoc__ = {}
__pdoc__["ateams.arithmetic.linearAlgebra"] = False
__pdoc__["ateams.arithmetic.linalg.KernelBasis"] = False
__pdoc__["ateams.arithmetic.linalg.SampleFromKernel"] = False
__pdoc__["ateams.arithmetic.cubicalComplex"] = False
__pdoc__["ateams.arithmetic.fastiteration"] = False
__pdoc__["ateams.arithmetic.reindexing"] = False

__pdoc__["ateams.arithmetic.common.MINT"] = "Global \"machine integer\" data type; this is equivalent to the C `int` data type."
__pdoc__["ateams.arithmetic.common.FINT"] = "Global \"finite field integer\" data type"
__pdoc__["ateams.arithmetic.common.SFINT"] = False

__all__ = [
	"Persistence", "MatrixReduction", "Kernel", "KernelSample"
]

