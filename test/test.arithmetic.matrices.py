
from ateams.arithmetic import SparseKernelBasis, KernelBasis, FINT
from galois import GF
import numpy as np
import sys

N = 20
p = 5
F = GF(p)
shapes = [(50,60), (60,50), (50,50)]

# Test over a series of random(ized) matrices of shape `shape`.
for shape in shapes:
	for _ in range(N):
		A = F.Random(shape)

		I = np.identity(shape[1], dtype=FINT)
		B = np.concatenate([A.T, I], axis=1, dtype=FINT)

		AUGMENT = shape[0]

		# Create arithmetic lookup tables.
		addition = np.zeros((p,p), dtype=FINT)
		for j in range(p): addition[:,j] = (np.arange(p, dtype=FINT)+j)%p

		subtraction = np.zeros((p,p), dtype=FINT)
		for j in range(p): subtraction[:,j] = (np.arange(p, dtype=FINT)-j)%p

		multiplication = np.zeros((p,p), dtype=FINT)
		for j in range(p): multiplication[:,j] = (np.arange(p, dtype=FINT)*j)%p

		negation = np.array([-q%p for q in range(0, p)], dtype=FINT)
		inverses = np.array([0] + list(GF(p).Range(1, p)**(-1))).astype(FINT)

		pivots = np.zeros(A.shape[1], dtype=FINT)
		sparse = bool(int(sys.argv[-2]))
		parallel = bool(int(sys.argv[-1]))
		empty = np.empty(shape=(0,0), dtype=FINT)
		
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
				assert np.array_equal(K, C)
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