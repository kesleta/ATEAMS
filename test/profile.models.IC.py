
import pstats, cProfile
import pyximport
import IC
import sys

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

M = IC.construct(L, sparse, parallel, minBlockSize, maxBlockSize, cores)

cProfile.runctx("IC.chain(M)", globals(), locals(), "Profile.prof")

fname = f"./profiles/IC/{L}{'.sparse' if sparse else ''}{'.parallel' if parallel else ''}.txt"
with open(fname, 'w') as stream:
    s = pstats.Stats('Profile.prof', stream=stream)
    s.strip_dirs().sort_stats("time").print_stats()

