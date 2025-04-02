
from ateams.arithmetic import SparseKernelBasis, KernelBasis
from galois import GF
import numpy as np
import sys

N = 20
p = 3
F = GF(p)
DTYPE = np.int64
shape = (50,60)

# Run a testing gamut: assert that our computations are correct on all possible
# 3x3 matrices over F_3.

# A = np.array([
# 	[0, 1, 1, 0, 1],
# 	[1, 2, 0, 1, 0],
# 	[0, 0, 0, 1, 2]
# ])
# A = F(A)
# shape = A.shape

# Test over a series of random(ized) matrices of shape `shape`.
for _ in range(N):
	A = F.Random(shape)

	I = np.identity(shape[1], dtype=DTYPE)
	B = np.concatenate([A.T, I], axis=1, dtype=DTYPE)

	AUGMENT = shape[0]

	# Create arithmetic lookup tables.
	addition = np.zeros((p,p), dtype=DTYPE)
	for j in range(p): addition[:,j] = (np.arange(p, dtype=DTYPE)+j)%p

	subtraction = np.zeros((p,p), dtype=DTYPE)
	for j in range(p): subtraction[:,j] = (np.arange(p, dtype=DTYPE)-j)%p

	multiplication = np.zeros((p,p), dtype=DTYPE)
	for j in range(p): multiplication[:,j] = (np.arange(p, dtype=DTYPE)*j)%p

	negation = np.array([-q%p for q in range(0, p)], dtype=DTYPE)
	inverses = np.array([0] + list(GF(p).Range(1, p)**(-1))).astype(DTYPE)

	pivots = np.zeros(A.shape[1], dtype=DTYPE)
	sparse = bool(int(sys.argv[-2]))
	parallel = bool(int(sys.argv[-1]))
	empty = np.empty(shape=(0,0), dtype=DTYPE)

	if sparse: C = np.asarray(SparseKernelBasis(pivots, empty, addition, subtraction, negation, multiplication, inverses, B, AUGMENT, parallel, 16, 32, 2))
	else: C = np.asarray(KernelBasis(pivots, empty, addition, subtraction, negation, multiplication, inverses, B, AUGMENT, parallel, 16, 32, 2))
		
	K = A.null_space()

	if K.shape[0] < 1 or K.shape[1] < 1:
		try:
			assert K.shape[0] == C.shape[0]
		except:
			print("###############################")
			print()
			print(C)
			print()
			print(K)
			print()
			print("###############################")
			sys.exit(1)
	else:
		try:
			assert np.array_equal(C, K)
		except:
			print("###############################")
			print()
			print(C)
			print()
			print(K)
			print()
			print("###############################")
			sys.exit(1)

sys.exit(0)