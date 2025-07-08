
import pstats, cProfile
import pyximport
import Persistence
import sys
import json


try:
    L = int(sys.argv[-6])
    sparse = bool(int(sys.argv[-5]))
    parallel = bool(int(sys.argv[-4]))
    minBlockSize = int(sys.argv[-3])
    maxBlockSize = int(sys.argv[-2])
    cores = int(sys.argv[-1])
except:
    L = 3
    sparse = False
    parallel = False
    minBlockSize = 32
    maxBlockSize = 64
    cores = 2
	
pyximport.install()

GROUND, TEST = Persistence.constructDefaults(L, sparse, parallel, minBlockSize, maxBlockSize, cores);

if sparse: PERSISTENCE = TEST
else: PERSISTENCE = GROUND

# Import test dataset.
with open(".testset.json") as r: TESTSET = json.load(r)
cProfile.runctx("Persistence.test(PERSISTENCE, TESTSET)", globals(), locals(), "Profile.prof")

fname = f"./profiles/persistence/{L}{'.sparse' if sparse else ''}{'.parallel' if parallel else ''}.txt"
with open(fname, 'w') as stream:
    s = pstats.Stats(pr, stream=stream)
    s.strip_dirs().sort_stats("time").print_stats()
