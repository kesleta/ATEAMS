
import warnings
warnings.simplefilter("ignore", UserWarning)

import pstats, cProfile
import pyximport
import IC
import sys


try:
    L = int(sys.argv[-7])
    sparse = bool(int(sys.argv[-6]))
    parallel = bool(int(sys.argv[-5]))
    minBlockSize = int(sys.argv[-4])
    maxBlockSize = int(sys.argv[-3])
    cores = int(sys.argv[-2])
    slurm = bool(int(sys.argv[-1]))
except:
    L = 3
    sparse = False
    parallel = False
    minBlockSize = 32
    maxBlockSize = 64
    cores = 2
    slurm = False

pyximport.install()

M = IC.construct(L, sparse, parallel, minBlockSize, maxBlockSize, cores)

TESTS = [L, sparse, parallel, minBlockSize, maxBlockSize, cores]
WIDTHS = [5, 10, 10, 5, 5, 5]
DESC = [str(int(thing)).ljust(width) for thing, width in zip(TESTS, WIDTHS)]
DESC = " ".join(DESC)

cProfile.runctx("IC.chain(M, DESC)", globals(), locals(), "Profile.prof")

fname = f"./profiles/InvadedCluster/{L}{'.sparse' if sparse else ''}{'.parallel' if parallel else ''}{'.slurm' if slurm else ''}.txt"
with open(fname, 'w') as stream:
    s = pstats.Stats('Profile.prof', stream=stream)
    s.strip_dirs().sort_stats("time").print_stats()

