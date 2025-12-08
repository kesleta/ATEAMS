import time
from ateams.complexes.Permutahedral import Permutahedral
from ateams.complexes.Cubical import Cubical
from ateams.models import Bernoulli
from ateams.models import SwendsenWang
from ateams.models import InvadedCluster
from ateams import Chain

start = time.time()

dimensions = [2,1]

# cmplx = Cubical().fromCorners(dimensions, periodic = False)
cmplx = Permutahedral().construct(dimensions, periodic = False)
# HP = Bernoulli(cmplx, field=2)
HP = SwendsenWang(cmplx, field=2)
# HP = InvadedCluster(cmplx, field=2)

print(cmplx.matrices.full.tolist())

# print(cmplx)
print(f"Seconds: {(time.time() - start):.2f}")

def dbg(x):
    print(x)
    return x

chain = dbg(Chain(HP, steps=10))

for x in chain:
    print(x)