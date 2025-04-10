
from ateams.arithmetic import Kernel, FINT, MatrixReduction
from galois import GF
import numpy as np
import sys
from tqdm import tqdm

N = 20
p = 3
F = GF(p)
shapes = [(50,60), (60,50), (50,50)]
# shapes = [(3,5), (5,3), (3,3)]
# shapes = [(9,12), (12,9), (12,12)]
# shapes = [(500,500)]

try:
	sparse = bool(int(sys.argv[-2]))
	parallel = bool(int(sys.argv[-1]))
except:
	sparse = True
	parallel = True

sp = "SPARSE" if sparse else "DENSE"
par = "PARALLEL" if parallel else "SERIAL"
DESC = (f"          {sp}, {par}").ljust(30)

# Test over a series of random(ized) matrices of shape `shape`.
with tqdm(total=len(shapes)*N, dynamic_ncols=True, desc=DESC) as bar:
	for shape in shapes:
		for _ in range(N):
			A = F.Random(shape)

			I = np.identity(shape[1], dtype=FINT)
			B = np.concatenate([A.T, I], axis=1, dtype=FINT)

			AUGMENT = shape[0]

			print()
			print()
			print("#################################")
			print("#################################")
			print()
			print(A)
			print()
			print(A.T)
			print()
			print("######## BEGIN REDUCTION")
			print()

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
			empty = np.empty(shape=(0,0), dtype=FINT)

			M = MatrixReduction(p, parallel, 16, 32, 2)
			
			# if sparse: C = np.asarray(SparseKernelBasisReduced(empty, M, AUGMENT))
			# else: C = np.asarray(SparseKernelBasis(pivots, empty, addition, subtraction, negation, multiplication, inverses, B, AUGMENT, parallel, 16, 32, 2))
			A = np.random.randint(0, p, size=shape, dtype=FINT)
			C = Kernel(M, A)
				
			K = F(A).null_space()

			print("######## END REDUCTION")

			if K.shape[0] < 1 or K.shape[1] < 1:
				flag = False;
				try:
					assert K.shape[0] == C.shape[0]
					print()
					print(C)
					print()
					print(K)
					print()
					print("#################################")
					print("#################################")
				except AssertionError:
					print()
					print(C)
					print()
					print(K)
					print()
					print("BAD")
					print("#################################")
					print("#################################")
					exit(1)
			else:
				try:
					assert np.array_equal(K, C)
					print()
					print(C)
					print()
					print(K)
					print()
					print("#################################")
					print("#################################")
				except AssertionError:
					print()
					print(C)
					print()
					print(K)
					print()
					print("BAD")
					print("#################################")
					print("#################################")
					exit(1)

			bar.update()

exit(0)